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

typedef struct HuffmanNode {
    HuffmanChar character;
    int frequency;
    bool isLeaf;
    struct HuffmanNode *left;
    struct HuffmanNode* right;
} HuffmanNode;

typedef struct PriorityQueue {
    long long int size;
    int capacity;
    HuffmanNode** list;
} PriorityQueue;

int THREAD_COUNT;

HuffmanChar* fin;
HuffmanEncoding* fout;

int inputFileSize;

int charFreq[CHAR_COUNT] = {0};
pthread_mutex_t charFreqMutex;
int** frequencyThread;

pthread_attr_t threadAttr;

typedef struct HuffmanCode {
    int encodingSize;
    HuffmanEncoding* code;
} HuffmanCode;

HuffmanCode huffmanCodeTable[CHAR_COUNT];

int *outputOffset, *outputLength;

int *tidList;

pthread_mutex_t* outputMutex;
pthread_mutex_t foutMutex;
pthread_cond_t* outputCond;

void *computeFreq(void *args);
void minHeapify(PriorityQueue *pq, int index);
HuffmanNode *extractMinPQ(PriorityQueue *pq);
int insertInPQ(PriorityQueue *pq, HuffmanNode *node);
void getHuffmanCodeLengths(HuffmanCode huffmanCodeTable[CHAR_COUNT], HuffmanNode *node, int len, int *max);
void getHuffmanCodes(HuffmanCode huffmanCodeTable[CHAR_COUNT], HuffmanNode *node, HuffmanEncoding *code, int len);
void writeBitsToBuffer(HuffmanEncoding *buf, int *offset, HuffmanEncoding code, int encodingSize);
void *encodeBatch(void *arg);
void cleanup(int outputFileSize, int inputFileSize, int outputFD, int inputFD);

void* computeFreq(void* args) {
    int tid = *((int*) args);
    int batchSize = (inputFileSize / THREAD_COUNT);
    if(tid == THREAD_COUNT - 1)
        batchSize = inputFileSize - (THREAD_COUNT - 1) * (inputFileSize / THREAD_COUNT);
    int startOffset = tid * (inputFileSize / THREAD_COUNT);
    for(int i = 0; i < batchSize; i++) 
        frequencyThread[tid][fin[startOffset + i]]++;
    pthread_mutex_lock(&charFreqMutex);
    for(int i = 0; i < CHAR_COUNT; i++){
        int temp = frequencyThread[tid][i];
        charFreq[i] += temp;
    }
    pthread_mutex_unlock(&charFreqMutex);
    pthread_exit(0);
}

void* encodeBatch(void* arg) {
    int tid = *((int*) arg);

    int len = 0;
    for(int i = 0; i < CHAR_COUNT; i++) {
        int tmpLen = frequencyThread[tid][i]*huffmanCodeTable[i].encodingSize;
        len += tmpLen;
    }
    
    outputLength[tid] = 0;
    __sync_fetch_and_add(&outputLength[tid], len);

    // calculate outputOffset, signal next thread
    int currOpOffset = 0;
    if(tid != 0) {
        while(__sync_fetch_and_add(&outputOffset[tid-1], 0) == -1 || __sync_fetch_and_add(&outputLength[tid-1], 0) == -1);
        currOpOffset = 0;
        __sync_fetch_and_add(&currOpOffset, outputOffset[tid - 1]);
        __sync_fetch_and_add(&currOpOffset, outputLength[tid - 1]);
        outputOffset[tid] = 0;
        __sync_fetch_and_add(&outputOffset[tid], currOpOffset);
    } else {
        currOpOffset = 0;
        __sync_fetch_and_add(&currOpOffset, outputOffset[tid]);
    }

    int first = currOpOffset >> 5;
    int last = (currOpOffset + len - 1) >> 5;

    if(tid != THREAD_COUNT - 1) {
        pthread_mutex_lock(&outputMutex[tid]);
        pthread_cond_signal(&outputCond[tid]);
        pthread_mutex_unlock(&outputMutex[tid]);
    }

    int startOffset = tid * (inputFileSize / THREAD_COUNT);
    

    int batchSize = (inputFileSize / THREAD_COUNT);
    if(tid == THREAD_COUNT - 1)
        batchSize += inputFileSize % THREAD_COUNT;


    for(int i = 0; i < batchSize; i++) {
        bool lock = false;
        int currFirst = currOpOffset >> 5;
        int currLast = (currOpOffset + huffmanCodeTable[fin[startOffset + i]].encodingSize - 1) >> 5;
        if(currFirst == first || currLast == last) lock = true;
        HuffmanCode hc = huffmanCodeTable[fin[startOffset + i]];
        int sz = hc.encodingSize;
        int cnt = ceil((double)sz / 32);
        for(int i = 0; i < cnt; i++) {
            int modEnc = hc.encodingSize%(1<<5);
            if(modEnc != 0 && i == cnt - 1) {
                if(lock) pthread_mutex_lock(&foutMutex);
                writeBitsToBuffer(fout, &currOpOffset, hc.code[i], modEnc);
                if(lock) pthread_mutex_unlock(&foutMutex);                    
            } else {
                if(lock) pthread_mutex_lock(&foutMutex);
                writeBitsToBuffer(fout, &currOpOffset, hc.code[i], 32);
                if(lock) pthread_mutex_unlock(&foutMutex);
            }
        }
    }
    free(frequencyThread[tid]);
    pthread_exit(0);
}

HuffmanNode* extractMinPQ(PriorityQueue* pq) {
    HuffmanNode* temp = pq->list[0];
    pq->list[0] = pq->list[pq->size - 1];
    pq->list[pq->size-1] = temp;
    pq->size--;
    minHeapify(pq, 0);
    return temp;
}

void minHeapify(PriorityQueue *pq, int index) {
    int min_index = index;
    int right_child = (index<<1) + 2;
    int left_child = right_child - 1;
    bool left_flag = false, right_flag = false;
    left_flag = (left_child >= 0 && left_child < pq->size);
    right_flag = (right_child >= 0 && right_child < pq->size);
    if(left_flag && pq->list[left_child]->frequency < pq->list[min_index]->frequency)
        min_index = left_child;
    if(right_flag && pq->list[right_child]->frequency < pq->list[min_index]->frequency)
        min_index = right_child;
    if(min_index != index) {
        HuffmanNode* tmp = pq->list[index];
        pq->list[index] = pq->list[min_index];
        pq->list[min_index] = tmp;
        minHeapify(pq, min_index);
    }
}


int insertInPQ(PriorityQueue* q, HuffmanNode* node) {
    q->size++;
    q->list[q->size - 1] = node;
    int i = q->size - 1;
    int par = (i+1)/2-1;
    while(i > 0 && q->list[par]->frequency > q->list[i]->frequency) {
        HuffmanNode* temp = q->list[i];
        q->list[i] = q->list[par];
        q->list[par] = temp;
        i = par;    //Propagate up the Huffman Tree
        par = (i+1)/2 - 1;
    }
    return 1;
}


void getHuffmanCodes(HuffmanCode huffmanCodeTable[CHAR_COUNT], HuffmanNode* node, HuffmanEncoding* code, int len) {
    if(node->isLeaf) {
        int cnt = ceil((double) len/32);
        for(int i = 0; i < cnt; i++) {
            huffmanCodeTable[node->character].code[i] = code[i];
        }
    }
    else {
        int offset = len >> 5;
        if(node->right) {
            code[offset]=code[offset]*2 + 1;
            getHuffmanCodes(huffmanCodeTable, node->right, code, len+1);
            code[offset]/=2;
        }
        if(node->left) {
            code[offset]*=2;
            getHuffmanCodes(huffmanCodeTable, node->left, code, len+1);
            code[offset]/=2;
        }
    }
    return;
}

void cleanup(int outputFileSize, int inputFileSize, int outputFD, int inputFD) {
    munmap(fout, outputFileSize);
    close(outputFD);
    munmap(fin, inputFileSize);
    close(inputFD);
}

void getHuffmanCodeLengths(HuffmanCode huffmanCodeTable[CHAR_COUNT], HuffmanNode* node, int len, int* max) {
    if(node->isLeaf) {
        if(len > *max) *max = len;
        int cnt = ceil((double) len/32);
        huffmanCodeTable[node->character].encodingSize = len;
        huffmanCodeTable[node->character].code = malloc(cnt*sizeof(HuffmanEncoding));

    } else {
        if(node->left) {
            getHuffmanCodeLengths(huffmanCodeTable, node->left, len+1, max);
        }
        if(node->right) {
            getHuffmanCodeLengths(huffmanCodeTable, node->right, len+1, max);
        }
    }
    return;
}

void writeBitsToBuffer(HuffmanEncoding* buf, int* offset, HuffmanEncoding code, int encodingSize) {
    int listOffset = *offset>>5;
    int extraOffset = *offset%(1<<5);
    int lim = (1<<5) - encodingSize - extraOffset;
    if(lim >= 0) {
        buf[listOffset] |= code << lim;
    } else {
        lim = -lim;
        buf[listOffset] |= code >> lim;
        buf[listOffset + 1] |= (code & ((1<<lim) - 1)) <<((1<<5)-lim);
    }
    *offset += encodingSize;
}

void frequencyUsingThreads(pthread_t *threads) {
    for(int i = 0; i < THREAD_COUNT; i++) {
        tidList[i] = i;
        frequencyThread[i] = malloc(CHAR_COUNT*sizeof(int));
        memset(frequencyThread[i], 0, CHAR_COUNT*sizeof(int));
        pthread_create(&threads[i], &threadAttr, computeFreq, &tidList[i]);
    }
    for(int i = 0; i < THREAD_COUNT; i++) {
        pthread_join(threads[i], NULL);
    }
}

void encodeUsingThreads(pthread_t *threads) {
    for(int i = 0; i < THREAD_COUNT; i++) {
        pthread_create(&threads[i], &threadAttr, encodeBatch, &tidList[i]);
    }
    for(int i = 0; i < THREAD_COUNT; i++) {
        pthread_join(threads[i], NULL);
    }
}

int main(int argc, char* argv[]) {
    double start, end;
    struct timeval timecheck;
    if(argc != 4) {
        fprintf(stderr, "Usage: %s <input text file> <output file name> <number of threads>\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    char* input = argv[1];
    char* output = argv[2];
    THREAD_COUNT = atoi(argv[3]);

    gettimeofday(&timecheck, NULL);
    start = timecheck.tv_sec + (double)timecheck.tv_usec / 1000000;

    // map input file
    int inputFD = open(input, O_RDONLY);
    inputFileSize = lseek(inputFD, 0, SEEK_END);
    lseek(inputFD, 0, SEEK_SET);
    fin = (HuffmanChar*) mmap(NULL, inputFileSize, PROT_READ, MAP_PRIVATE, inputFD, 0);

    // measure char frequency for ASCII set
    pthread_mutex_init(&charFreqMutex, NULL);
    frequencyThread = malloc(sizeof(int*) * THREAD_COUNT);

    outputLength = (int*) malloc(THREAD_COUNT*sizeof(int));

    // initialize thread variables
    pthread_t threads[THREAD_COUNT];
    pthread_attr_init(&threadAttr);

    tidList = malloc(sizeof(int) * THREAD_COUNT);
    
    frequencyUsingThreads(threads);

    // construct huffman tree
    PriorityQueue* pq = malloc(sizeof(PriorityQueue));
    pq->capacity = 256;
    for(int i = 0; i < CHAR_COUNT; i++) 
        if(charFreq[i] > 0)
            pq->size++;
    
    pq->list = malloc(pq->size * sizeof(HuffmanNode*));

    int j = 0;
    for(int i = 0; i < CHAR_COUNT; i++) {
        if(charFreq[i] > 0) {
            pq->list[j] = malloc(sizeof(HuffmanNode));
            pq->list[j]->left = NULL;
            pq->list[j]->right = NULL;
            pq->list[j]->isLeaf = true;
            pq->list[j]->character = i;
            pq->list[j]->frequency = charFreq[i];
            j++;
        }
    }
    int num = pq->size;
    for(int i = num / 2 - 1; i >= 0; i--) {
        minHeapify(pq, i);
    }
    int n = pq->size;
    for(int i = 0; i < n - 1; i++) {
        HuffmanNode* left = extractMinPQ(pq);
        HuffmanNode* right = extractMinPQ(pq);
        HuffmanNode* node = malloc(sizeof(HuffmanNode));
        node->left = left;
        node->right = right;
        node->isLeaf = false;
        node->character = 0;
        node->frequency = left->frequency + right->frequency;
        insertInPQ(pq, node);
    }


    outputCond = (pthread_cond_t*) malloc(THREAD_COUNT*sizeof(pthread_cond_t));

    HuffmanNode* root = extractMinPQ(pq);

    for(int i = 0; i < CHAR_COUNT; i++) {
        huffmanCodeTable[i].encodingSize = -1;
    }

    int maxlen = 0;
    getHuffmanCodeLengths(huffmanCodeTable, root, 0, &maxlen);
    HuffmanEncoding* codeval = malloc(sizeof(HuffmanEncoding) * (ceil((double) maxlen / 32)));
    getHuffmanCodes(huffmanCodeTable, root, codeval, 0);

    int outputFileSize = 72 + 32 * THREAD_COUNT;
    for(int i = 0; i < CHAR_COUNT; i++) {
        if(huffmanCodeTable[i].encodingSize > 0) {
            int sz = 16+huffmanCodeTable[i].encodingSize;
            outputFileSize += sz;
        }
        outputFileSize += charFreq[i] * huffmanCodeTable[i].encodingSize;
    }
    outputFileSize = ceil((double)outputFileSize/32) * sizeof(HuffmanEncoding);
    int outputFD = open(output, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH);
    ftruncate(outputFD, outputFileSize);


    outputMutex = (pthread_mutex_t*) malloc(THREAD_COUNT*sizeof(pthread_mutex_t));

    fout = (HuffmanEncoding*) mmap(NULL, outputFileSize, PROT_READ | PROT_WRITE, MAP_SHARED, outputFD, 0);
    memset(fout, 0, outputFileSize);

    int bitOffset = 0;
    writeBitsToBuffer(fout, &bitOffset, (HuffmanEncoding)inputFileSize, 32);
    writeBitsToBuffer(fout, &bitOffset, (HuffmanEncoding)THREAD_COUNT, 32);
    int hdr = 0;
    while(hdr < THREAD_COUNT) {
        writeBitsToBuffer(fout, &bitOffset, (HuffmanEncoding)0, 32);
        hdr++;
    }

    HuffmanEncoding counter = 0;
    for(int i = 0; i < CHAR_COUNT; i++)
        if(huffmanCodeTable[i].encodingSize > 0)
            counter++;
    writeBitsToBuffer(fout, &bitOffset, counter, 8);
    for(int i = 0; i < CHAR_COUNT; i++) {
        if(huffmanCodeTable[i].encodingSize > 0) {
            writeBitsToBuffer(fout, &bitOffset, (HuffmanEncoding)i, 8);
            writeBitsToBuffer(fout, &bitOffset, (HuffmanEncoding)huffmanCodeTable[i].encodingSize, 8);
            HuffmanCode hc = huffmanCodeTable[i];
            int cnt = ceil((double)hc.encodingSize/(1<<5));
            for(int i = 0; i < cnt; i++) {
                int modEnc = hc.encodingSize%(1<<5);
                if(modEnc != 0 && i == cnt - 1)
                    writeBitsToBuffer(fout, &bitOffset, hc.code[i], modEnc);
                else writeBitsToBuffer(fout, &bitOffset, hc.code[i], 32);
            }
        }
    }
    
    outputOffset = (int*) malloc(THREAD_COUNT*sizeof(int));
    pthread_mutex_init(&foutMutex, NULL);
    for(int i = 0; i < THREAD_COUNT; i++) {
        outputLength[i] = outputOffset[i] = -1;
        pthread_mutex_init(&outputMutex[i], NULL);
        pthread_cond_init(&outputCond[i], NULL);
    }    
    outputOffset[0] = bitOffset;
    encodeUsingThreads(threads);

    for(int i = 2; i<=THREAD_COUNT+1; i++) {
        fout[i] = (HuffmanEncoding)outputOffset[i-2];
    }

    cleanup(outputFileSize, inputFileSize, outputFD, inputFD);

    gettimeofday(&timecheck, NULL);
    end = timecheck.tv_sec + (double) timecheck.tv_usec / 1000000;
    printf("Encoding: %lf sec\n", end - start);
}