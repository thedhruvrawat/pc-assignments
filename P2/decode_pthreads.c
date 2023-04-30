#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include <sys/time.h>
#include <stdbool.h>
#include <stdatomic.h>
#include <sys/stat.h>


#define CHAR_COUNT 256

typedef uint32_t HuffmanEncoding;
typedef uint8_t HuffmanChar;

typedef struct HuffmanTreeNode {
    HuffmanChar character;
    bool isLeaf;
    struct HuffmanTreeNode *left;
    struct HuffmanTreeNode* right;
} HuffmanTreeNode;

int THREAD_COUNT;

HuffmanEncoding* fin;
HuffmanChar* fout;

int outputFileSize;

typedef struct HuffmanCode {
    int encodingSize;
    HuffmanEncoding* code;
} HuffmanCode;

int* inputOffset;

HuffmanTreeNode* root;

HuffmanEncoding readBitsFromBuffer(HuffmanEncoding *buf, int *offset, int encodingSize);
void *decodeBatch(void *args);
void cleanup(int outputFileSize, int inputFileSize, int outputFD, int inputFD);
int readInput(char *input, int *inputFileSize);
int writeOutput(char *output, int *outputFD);

void cleanup(int outputFileSize, int inputFileSize, int outputFD, int inputFD)
{
    munmap(fout, outputFileSize);
    close(outputFD);
    munmap(fin, inputFileSize);
    close(inputFD);
}

HuffmanEncoding readBitsFromBuffer(HuffmanEncoding* buf, int* offset, int encodingSize) {
    HuffmanEncoding res = 0;
    int off = *offset>>5;
    int rem = (1<<5) - (*offset%(1<<5)) - 1;
    *offset += encodingSize;
    HuffmanEncoding cover = 1 << rem;
    for(int i = 0; i < encodingSize; i++) {
        res *= 2;
        if (buf[off] & cover) res++;
        cover/=2;
        if(cover == 0) {
            off++;
            cover = 1 << ((1<<5)-1);
        }
    }
    return res;
}

void* decodeBatch(void* args) {
    int tid = *((int*) args);
    int batchSize = (outputFileSize / THREAD_COUNT);
    if(tid == THREAD_COUNT - 1)
        batchSize += outputFileSize % THREAD_COUNT;
    int off = inputOffset[tid];
    for(int i = 0; i < batchSize; i++) {
        int listOffset = off>>5;
        int extraOffset = off%(1<<5);
        HuffmanEncoding cover = 1 << ((1<<5) - extraOffset - 1);
        HuffmanChar hc = 0;
        int len = 0;
        HuffmanTreeNode* node = root;
        for (;;) {
            bool bit = (fin[listOffset] & cover);
            cover/=2;
            if(!bit) node = node->left;
            else node = node->right;
            len++;
            if(node == NULL) {
                fprintf(stderr, "%s\n", "Error during Decoding");
            }
            if(node->isLeaf) {
                hc = node->character;
                break;
            }
            if(cover == 0) {
                listOffset++;
                cover = 1<<((1<<5)-1);
            }
        }
        off += len;
        fout[tid * (outputFileSize / THREAD_COUNT) + i] = hc;
    }
    pthread_exit(0);
}

int readInput(char *input, int *inputFileSize) {
    int inputFD = open(input, O_RDONLY);
    int filesize = lseek(inputFD, 0, SEEK_END);
    *inputFileSize = filesize;
    lseek(inputFD, 0, SEEK_SET);
    fin = (HuffmanEncoding*) mmap(NULL, filesize, PROT_READ, MAP_PRIVATE, inputFD, 0);
    return inputFD;
}

int writeOutput(char *output, int *outputFD) {
    int off = 0;
    outputFileSize = readBitsFromBuffer(fin, &off, 32);
    int fd = open(output, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH);
    *outputFD = fd;
    ftruncate(fd, outputFileSize);
    fout = (HuffmanChar*) mmap(NULL, outputFileSize, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    return off;
}

void decodeUsingThreads(int *tidList) {
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_t threads[THREAD_COUNT];
    for(int i = 0; i < THREAD_COUNT; i++) {
        pthread_create(&threads[i], &attr, decodeBatch, &tidList[i]);
    }
    for(int i = 0; i < THREAD_COUNT; i++) {
        pthread_join(threads[i], NULL);
    }
}

int main(int argc, char* argv[]) {
    double start, end;
    struct timeval timecheck;

    if(argc != 3) {
        fprintf(stderr, "Usage: %s <encoded text file> <output file name>\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    char* input = argv[1];
    char* output = argv[2];

    gettimeofday(&timecheck, NULL);
    start = timecheck.tv_sec + (double)timecheck.tv_usec/1000000;

    int inputFileSize = 0;

    int inputFD = readInput(input, &inputFileSize);
    int outputFD = -1;
    int bitOffset = writeOutput(output, &outputFD);

    HuffmanCode codeTable[CHAR_COUNT];
    memset(fout, 0, outputFileSize);

    for(int i = 0; i < CHAR_COUNT; i++)
        codeTable[i].encodingSize = -1;

    THREAD_COUNT = readBitsFromBuffer(fin, &bitOffset, 32);

    int *tidList = (int*) malloc(sizeof(int) * THREAD_COUNT);
    for(int i = 0; i < THREAD_COUNT; i++) {
        tidList[i] = i;
    }

    inputOffset = (int*) malloc(sizeof(int) * THREAD_COUNT);
    for(int i = 0; i < THREAD_COUNT; i++) {
        inputOffset[i] = readBitsFromBuffer(fin, &bitOffset, 32);
    }

    HuffmanTreeNode* node = malloc(sizeof(HuffmanTreeNode));
    node->character = 0;
    node->left = NULL;
    node->right = NULL;
    node->isLeaf = true;
    root = node;


    HuffmanEncoding char_cnt = readBitsFromBuffer(fin, &bitOffset, 8);
    for(int i = 0; i < char_cnt; i++) {
        HuffmanEncoding hc = readBitsFromBuffer(fin, &bitOffset, 8);
        int sz = (int)readBitsFromBuffer(fin, &bitOffset, 8);
        codeTable[hc].encodingSize = sz;
        int cnt = ceil((double) sz/(1<<5));
        int modEnc = sz % (1<<5);
        codeTable[hc].code = malloc(cnt*sizeof(HuffmanEncoding));
        for(int i = 0; i < cnt; i++) {
            if(modEnc != 0 && i == cnt - 1) {
                codeTable[hc].code[i] = readBitsFromBuffer(fin, &bitOffset, modEnc);
            } else {
                codeTable[hc].code[i] = readBitsFromBuffer(fin, &bitOffset, 32);
            }
        }
    }   
    
    for (int i = 0; i < CHAR_COUNT; i++) {
        int sz = codeTable[i].encodingSize;
        if (sz > 0) {
            HuffmanTreeNode *temp = root;
            int cnt = ceil((double)sz/(1<<5));
            for(int j = 0; j < cnt; j++) {
                int modEnc = sz % (1<<5);
                modEnc--;
                HuffmanEncoding cover;
                if(modEnc+1 != 0 && j == cnt - 1)
                    cover = 1 << modEnc;
                else {
                    cover = 1 << ((1<<5) - 1);
                }
                while(cover != 0) {
                    if(codeTable[i].code[j] & cover) {
                        if(temp->right == NULL) {
                            HuffmanTreeNode* node = malloc(sizeof(HuffmanTreeNode));
                            node->character = 0;
                            node->left = NULL;
                            node->right = NULL;
                            node->isLeaf = true;
                            temp->right = node;
                        }
                        temp->isLeaf = false;
                        temp = temp->right;
                    } else {
                        if(temp->left == NULL) {
                            HuffmanTreeNode* node = malloc(sizeof(HuffmanTreeNode));
                            node->character = 0;
                            node->left = NULL;
                            node->right = NULL;
                            node->isLeaf = true;
                            temp->left = node;
                        }
                        temp->isLeaf = false;
                        temp = temp->left;
                    }
                    cover/=2;
                }
            }
            temp->character = i;
        }
    }

    decodeUsingThreads(tidList);

    gettimeofday(&timecheck, NULL);
    end = timecheck.tv_sec + (double) timecheck.tv_usec / 1000000;
    cleanup(outputFileSize, inputFileSize, outputFD, inputFD);
    printf("Decoding: %lf sec\n", end - start);
}