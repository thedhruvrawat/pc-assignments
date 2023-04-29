#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <omp.h> 
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#include <unistd.h>
#include <time.h>

#define MAX_HEAP_SIZE 256

#define MAX_DECODED_BUFFER_SIZE 4096
#define MAX_ENCODED_BUFFER_SIZE 32*MAX_DECODED_BUFFER_SIZE
#define PATH_MAX 4096

typedef struct Dictionary {
    unsigned long long frequencies[MAX_HEAP_SIZE];
    int size;
} Dictionary;

typedef struct Node {
    int value;  // 0-255 for data, -1 for placeholder
    unsigned long long priority;
    struct Node* left;
    struct Node* right;
} Node;

typedef struct PriorityQ {
    Node** minHeap;
    int size;
    int capacity;
} PriorityQ;

Dictionary* createDictionary(int size);
void freeDictionary(Dictionary* d);
void mergeDictionaries(Dictionary* d, Dictionary** dictionaries, int num_threads);

Node* createNode(int value, unsigned long long priority, Node* left, Node* right);
void freeNode(Node* node);
PriorityQ* createPriorityQ();
void freePriorityQ(PriorityQ* pq);
bool swapMinHeapElements(Node** minHeap, int index_a, int index_b);
bool pushPriorityQ(PriorityQ* pq, Node* node);
bool popPriorityQ(PriorityQ* pq, Node** node);

Node* getHuffmanTree(Dictionary* dict);
void getHuffmanAlphabet(Node* node, int level, char currentCode[], char* huffmanAlphabet[]);
void freeHuffmanAlphabet(char* huffmanAlphabet[]);
bool fwriteBitInBuffer(char bit, char buffer[], int* nbits, unsigned long long* nbytes);

bool chunkEncoder(unsigned char* inputChunk, unsigned char* outputChunk, char* huffmanAlphabet[], unsigned long long inputBufferSize, unsigned long long* outputChunkSize);
bool fileEncoderBarrier(FILE* inputFile, FILE* outputFile, char* huffmanAlphabet[], unsigned long long* outputFileSize, unsigned long long inputChunkSizes[], unsigned long long outputChunkSizes[], int num_threads);
bool fileEncoderLocks(FILE* inputFile, FILE* outputFile, char* huffmanAlphabet[], unsigned long long* outputFileSize, unsigned long long inputChunkSizes[], unsigned long long outputChunkSizes[], int num_threads);
bool fileEncoderFull(char* inputFileName, char* outputFileName, int num_threads, int mode, int rank);

void countChunk(unsigned char* chunk, unsigned long long size, Dictionary* d);
unsigned long long parallel_get_frequencies(FILE* file, Dictionary* d, int num_threads);

Dictionary* createDictionary(int size) {
    Dictionary* d = (Dictionary*)malloc(sizeof(Dictionary));
    for (int i = 0; i < size; i++) {
        d->frequencies[i] = 0;
    }
    d->size = size;
    return d;
}

void freeDictionary(Dictionary* d) {
    free(d);
    return;
}

void mergeDictionaries(Dictionary* d, Dictionary** dictionaries, int num_threads) {
    for (int i = 0; i < num_threads; i++) {
        for (int j = 0; j < d->size; j++) {
            d->frequencies[j] += dictionaries[i]->frequencies[j];
        }
    }
}

Node* createNode(int value, unsigned long long priority, Node* left, Node* right) {
    Node* node = (Node*)malloc(sizeof(Node));
    node->left = left;
    node->right = right;
    node->value = value;
    node->priority = priority;
    return node;
}

void freeNode(Node* node){
    if (node->left != NULL){
        freeNode(node->left);
        node->left = NULL;
    }
    if (node->right != NULL){
        freeNode(node->right);
        node->right = NULL;
    }
    free(node);
}

PriorityQ* createPriorityQ() {
    PriorityQ* pq = (PriorityQ*)malloc(sizeof(PriorityQ));
    pq->minHeap = (Node**)malloc(sizeof(Node*) * MAX_HEAP_SIZE);
    pq->size = 0;
    pq->capacity = MAX_HEAP_SIZE;
    // pq->minHeap is an array of Node* stored in the stack
    return pq;
}

void freePriorityQ(PriorityQ* pq) {
    for (int i = 0; i < pq->size; i++) {
        free(pq->minHeap[i]);
        pq->minHeap[i] = NULL;
    }
    free(pq->minHeap);
    free(pq);
    return;
}

bool swapMinHeapElements(Node** minHeap, int index_a, int index_b) {
    Node* temp = *(minHeap + index_a);
    *(minHeap + index_a) = *(minHeap + index_b);
    *(minHeap + index_b) = temp;
    return true;
}

bool pushPriorityQ(PriorityQ* pq, Node* node) {
    if (pq->size == pq->capacity) {
        return false;
    }
    pq->minHeap[pq->size] = node;
    pq->size++;
    int i = pq->size - 1;
    while (i > 0 && pq->minHeap[i]->priority < pq->minHeap[(i - 1) / 2]->priority) {
        swapMinHeapElements(pq->minHeap, i, (i - 1) / 2);
        i = (i - 1) / 2;
    }
    return true;
}

bool popPriorityQ(PriorityQ* pq, Node** node) {
    if (pq->size == 0) {
        return false;
    }
    // pop the first element, which is the minimum
    *node = pq->minHeap[0];
    pq->minHeap[0] = pq->minHeap[pq->size - 1];
    pq->size--;
    int i = 0;
    while (i < pq->size) {
        Node* current = pq->minHeap[i];
        int leftChildIndex = 2 * i + 1;
        int rightChildIndex = 2 * i + 2;
        // if all children
        if (leftChildIndex < pq->size && leftChildIndex < pq->size) {
            Node* leftChild = pq->minHeap[leftChildIndex];
            Node* rightChild = pq->minHeap[rightChildIndex];
            if (current->priority > leftChild->priority || current->priority > rightChild->priority) {
                if (leftChild->priority < rightChild->priority) {
                    swapMinHeapElements(pq->minHeap, i, leftChildIndex);
                    i = leftChildIndex;
                }
                else {
                    swapMinHeapElements(pq->minHeap, i, rightChildIndex);
                    i = rightChildIndex;
                }
            }
            else {
                break;
            }
        }
        // if no right child
        else if (2 * i + 1 < pq->size) {
            Node* leftChild = pq->minHeap[leftChildIndex];
            if (current->priority < leftChild->priority) {
                swapMinHeapElements(pq->minHeap, i, leftChildIndex);
                i = leftChildIndex;
            }
            else {
                break;
            }
        }
        // if no children
        else {
            break;
        }
    }
    return true;
}

Node* getHuffmanTree(Dictionary* dict) {
    PriorityQ* pq = createPriorityQ(MAX_HEAP_SIZE);
    for (int i = 0; i < dict->size; i++) {
        Node* node = createNode(i, dict->frequencies[i], NULL, NULL);
        pushPriorityQ(pq, node);
    }
    for (int i = 0; i < dict->size - 1; i++) {
        Node* child1;
        Node* child2;
        popPriorityQ(pq, &child1);
        popPriorityQ(pq, &child2);
        Node* newNode = createNode(-1, child1->priority + child2->priority, child1, child2);
        pushPriorityQ(pq, newNode);
    }
    Node* huffmanTree;
    popPriorityQ(pq, &huffmanTree);
    freePriorityQ(pq);
    pq = NULL;
    return huffmanTree;
}


void getHuffmanAlphabet(Node* node, int level, char currentCode[], char* huffmanAlphabet[]) {

    if ((node->left == NULL) && (node->right == NULL)) {
        currentCode[level] = 0;  // end of string
        huffmanAlphabet[node->value] = strdup(currentCode);
    }
    else {
        // left -> 0
        currentCode[level] = '0';
        getHuffmanAlphabet(node->left, level + 1, currentCode, huffmanAlphabet);
        // right -> 1
        currentCode[level] = '1';
        getHuffmanAlphabet(node->right, level + 1, currentCode, huffmanAlphabet);
    }
}

void freeHuffmanAlphabet(char* huffmanAlphabet[]) {
    for (int i = 0; i < MAX_HEAP_SIZE; i++) {
        free(huffmanAlphabet[i]);
        huffmanAlphabet[i] = NULL;
    }
}

bool fwriteBitInBuffer(char bit, char buffer[], int* nbits, unsigned long long* nbytes) {
    buffer[*nbytes] <<= 1;
    if (bit == '1') buffer[*nbytes] |= 1;
    (*nbits)++;
    if (*nbits == 8) {
        *nbits = 0;
        (*nbytes)++;
    }
    if (*nbytes == MAX_ENCODED_BUFFER_SIZE) {
        return false;
    }
    return true;
}

void countChunk(unsigned char* chunk, unsigned long long size, Dictionary* d) {
    for (unsigned long long i = 0; i < size; i++) {
        int value = chunk[i];
        #pragma omp atomic update
        d->frequencies[value]++;
    }
}

unsigned long long parallel_get_frequencies(FILE* file, Dictionary* d, int num_threads) {
    Dictionary** dictionaries = (Dictionary**)malloc(num_threads * sizeof(Dictionary*));
    for (int i = 0; i < num_threads; i++) {
        dictionaries[i] = createDictionary(MAX_HEAP_SIZE);
    }
    fseek(file, 0, SEEK_END);       // seek to end of file
    unsigned long long file_size = ftell(file);    // get current file pointer
    fseek(file, 0, SEEK_SET);       // seek to end of file
    // chunks
    unsigned char* chunk = (unsigned char*)malloc(sizeof(unsigned char) * num_threads * MAX_DECODED_BUFFER_SIZE);

    unsigned long long chunk_size = MAX_DECODED_BUFFER_SIZE;
    unsigned long long chunk_count = 0;
    unsigned long long* read = malloc(sizeof(unsigned long long) * num_threads);

    chunk_count = file_size / chunk_size;
    if (file_size % chunk_size != 0) {
        chunk_count++;
    }
    unsigned long long chunk_iterations = chunk_count / num_threads;
    if (chunk_count % num_threads != 0) {
        chunk_iterations++;
    }
    #pragma omp parallel
    for (unsigned long long i = 0; i < chunk_iterations; i++) {
        int thread_ID = 0;

        thread_ID = omp_get_thread_num();
        // master sequential read of num_threads chunks
        #pragma omp single
        {
            for (int j = 0;j < num_threads;j++) {
                read[j] = fread(chunk + j * MAX_DECODED_BUFFER_SIZE, sizeof(unsigned char), chunk_size, file);
            }
        }
        //  count the chunks
        if (read[thread_ID] > 0)
            countChunk(chunk + thread_ID * MAX_DECODED_BUFFER_SIZE, read[thread_ID], dictionaries[thread_ID]);

    }
    // merge dictionaries
    mergeDictionaries(d, dictionaries, num_threads);
    for (int i = 0; i < num_threads; i++) {
        freeDictionary(dictionaries[i]);
    }
    free(dictionaries);
    free(chunk);
    free(read);
    return file_size;
}

unsigned char getCharFromHuffmanEncodedBitStream(unsigned char buffer[], unsigned long long* nbytes, int* nbits, Node* node) {
    if (node->left == NULL && node->right == NULL) {
        return node->value;
    }

    bool isBitSet = buffer[*nbytes] & (1 << (7 - *nbits));

    (*nbits)++;
    if (*nbits == 8) {
        *nbits = 0;
        (*nbytes)++;
    }


    if (!isBitSet && node->left != NULL) {
        return getCharFromHuffmanEncodedBitStream(buffer, nbytes, nbits, node->left);
    }
    else if (node->right != NULL) {
        return getCharFromHuffmanEncodedBitStream(buffer, nbytes, nbits, node->right);
    }
    return '\0';
}

bool getFrequenciesFromEncodedFile(FILE* inputFile, Dictionary* dict) {
    fseek(inputFile, -(sizeof(unsigned long long) * (MAX_HEAP_SIZE + 1)), SEEK_END);

    int numFrequencies = fread(dict->frequencies, sizeof(unsigned long long), MAX_HEAP_SIZE, inputFile);

    return numFrequencies == MAX_HEAP_SIZE;
}

unsigned long long* getChunkOffsetsFromEncodedFile(FILE* inputFile, unsigned long long* numChunks) {
    fseek(inputFile, -sizeof(unsigned long long), SEEK_END);

    fread(numChunks, sizeof(unsigned long long), 1, inputFile);

    unsigned long long* chunkOffsets = (unsigned long long*)malloc(sizeof(unsigned long long) * (*numChunks + 1));
    fseek(inputFile, -(sizeof(unsigned long long) * (MAX_HEAP_SIZE + 2 * *numChunks + 1 + 1)), SEEK_END);
    fread(chunkOffsets, sizeof(unsigned long long), *numChunks + 1, inputFile);

    return chunkOffsets;
}

unsigned long long* getOriginalChunkSizesFromEncodedFile(FILE* inputFile, unsigned long long* numChunks) {
    fseek(inputFile, -sizeof(unsigned long long), SEEK_END);

    fread(numChunks, sizeof(unsigned long long), 1, inputFile);

    unsigned long long* chunkOffsets = (unsigned long long*)malloc(sizeof(unsigned long long) * (*numChunks));
    fseek(inputFile, -(sizeof(unsigned long long) * (MAX_HEAP_SIZE + *numChunks + 1)), SEEK_END);
    fread(chunkOffsets, sizeof(unsigned long long), *numChunks, inputFile);

    return chunkOffsets;
}

bool chunkDecoder(unsigned char inputChunk[], unsigned char outputChunk[], Node* huffmanTree, unsigned long long inputChunkSize, unsigned long long* outputChunkSize) {
    unsigned long long nbytes = 0;
    int nbits = 0;
    unsigned long long outputCharCounter = 0;
    bool isDecodingSuccessful = true;

    // decode the chunk
    while (outputCharCounter < inputChunkSize) {
        outputChunk[outputCharCounter] = getCharFromHuffmanEncodedBitStream(inputChunk, &nbytes, &nbits, huffmanTree);

        outputCharCounter++;
    }
    *outputChunkSize = outputCharCounter;
    return isDecodingSuccessful;
}

bool fileDecoderBarrier(FILE* inputFile, FILE* outputFile, Node* huffmanTree, unsigned long long inputChunkOffsets[], unsigned long long inputChunkSizes[], unsigned long long numOfChunks, int num_threads) {
    // chunks
    unsigned char* inputChunk = (unsigned char*)malloc(sizeof(unsigned char) * num_threads * MAX_ENCODED_BUFFER_SIZE);
    // +1 to avoid overflow
    unsigned char* outputChunk = (unsigned char*)malloc(sizeof(unsigned char) * num_threads * (MAX_DECODED_BUFFER_SIZE));
    bool isDecodingSuccessful = true;
    unsigned long long* inputBufferChunkSizes = (unsigned long long*)malloc(sizeof(unsigned long long) * num_threads);
    unsigned long long* outputBufferChunkSizes = (unsigned long long*)malloc(sizeof(unsigned long long) * num_threads);
    unsigned long long* decodedOutputBufferChunkSizes = (unsigned long long*)malloc(sizeof(unsigned long long) * num_threads);
    unsigned long long* inputBufferChunkOffsets = (unsigned long long*)malloc(sizeof(unsigned long long) * num_threads);
    unsigned long long chunkIterations = numOfChunks / num_threads;
    if (numOfChunks % num_threads != 0) {
        chunkIterations++;
    }
#pragma omp parallel
    for (unsigned long long i = 0; i < chunkIterations; i++) {
        int thread_ID = 0;
        thread_ID = omp_get_thread_num();
        // single master thread as reader
        #pragma omp single
        {
            unsigned long long readSize = 0;
            for (int j = 0;j < num_threads;j++) {
                if (i * num_threads + j < numOfChunks) {
                    inputBufferChunkSizes[j] = inputChunkOffsets[j + i * num_threads + 1] - inputChunkOffsets[j + i * num_threads];
                    outputBufferChunkSizes[j] = inputChunkSizes[j + i * num_threads];
                    inputBufferChunkOffsets[j] = inputChunkOffsets[j + i * num_threads];
                }
                else {
                    inputBufferChunkSizes[j] = 0;
                    outputBufferChunkSizes[j] = 0;
                    inputBufferChunkOffsets[j] = 0;
                }
                if (inputBufferChunkSizes[j] > 0) {
                    fseek(inputFile, inputBufferChunkOffsets[j], SEEK_SET);
                    readSize = fread(inputChunk + j * MAX_ENCODED_BUFFER_SIZE, sizeof(unsigned char), inputBufferChunkSizes[j], inputFile);
                    if (readSize != inputBufferChunkSizes[j]) {
                        isDecodingSuccessful = false;
                    }
                }
            }
        }
        if (outputBufferChunkSizes[thread_ID] > 0) {
            chunkDecoder(inputChunk + thread_ID * MAX_ENCODED_BUFFER_SIZE, outputChunk + thread_ID * (MAX_DECODED_BUFFER_SIZE), huffmanTree, outputBufferChunkSizes[thread_ID], &decodedOutputBufferChunkSizes[thread_ID]);
        }
#pragma omp barrier
#pragma omp single
        {
            // writer
            for (int j = 0;j < num_threads;j++) {
                if (outputBufferChunkSizes[j] > 0) {
                    fwrite(outputChunk + j * (MAX_DECODED_BUFFER_SIZE), sizeof(unsigned char), decodedOutputBufferChunkSizes[j], outputFile);
                }
            }
        }
    }
    free(inputChunk);
    free(outputChunk);
    free(inputBufferChunkSizes);
    free(outputBufferChunkSizes);
    free(decodedOutputBufferChunkSizes);
    free(inputBufferChunkOffsets);

    return isDecodingSuccessful;
}

bool fileDecoderLocks(FILE* inputFile, FILE* outputFile, Node* huffmanTree, unsigned long long inputChunkOffsets[], unsigned long long inputChunkSizes[], unsigned long long numOfChunks, int num_threads) {
    unsigned char* inputChunk = (unsigned char*)malloc(sizeof(unsigned char) * num_threads * MAX_ENCODED_BUFFER_SIZE);
    // +1 to avoid overflow
    unsigned char* outputChunk = (unsigned char*)malloc(sizeof(unsigned char) * num_threads * (MAX_DECODED_BUFFER_SIZE));
    bool isDecodingSuccessful = true;
    unsigned long long* inputBufferChunkSizes = (unsigned long long*)malloc(sizeof(unsigned long long) * num_threads);
    unsigned long long* outputBufferChunkSizes = (unsigned long long*)malloc(sizeof(unsigned long long) * num_threads);
    unsigned long long* decodedOutputBufferChunkSizes = (unsigned long long*)malloc(sizeof(unsigned long long) * num_threads);
    unsigned long long* inputBufferChunkOffsets = (unsigned long long*)malloc(sizeof(unsigned long long) * num_threads);
    unsigned long long chunkIterations = numOfChunks / num_threads;
    if (numOfChunks % num_threads != 0) {
        chunkIterations++;
    }
    omp_lock_t* readlock = malloc(sizeof(omp_lock_t) * num_threads);
    omp_lock_t* processlock = malloc(sizeof(omp_lock_t) * num_threads);
    omp_lock_t* writelock = malloc(sizeof(omp_lock_t) * num_threads);
    for (int i = 0;i < num_threads;i++) {
        omp_lock_t r_lock, p_lock, w_lock;
        readlock[i] = r_lock;
        processlock[i] = p_lock;
        writelock[i] = w_lock;
    }
    int current_chunk = 0;
    for (int j = 0;j < num_threads;j++) {
        omp_init_lock(&readlock[j]);
        omp_init_lock(&processlock[j]);
        omp_init_lock(&writelock[j]);
        omp_unset_lock(&readlock[j]);
        omp_set_lock(&processlock[j]);
        omp_set_lock(&writelock[j]);
    }
    omp_set_dynamic(0);
    omp_set_num_threads(num_threads + 2);
    #pragma omp parallel
    for (unsigned long long i = 0; i < chunkIterations; i++) {
        int thread_ID = omp_get_thread_num();
        // single master thread as reader
        if (thread_ID == 0) {
            unsigned long long readSize = 0;
            for (unsigned long long j = 0;j < num_threads;j++) {
                omp_set_lock(&readlock[j]);
                if (i * num_threads + j < numOfChunks) {
                    inputBufferChunkSizes[j] = inputChunkOffsets[j + i * num_threads + 1] - inputChunkOffsets[j + i * num_threads];
                    outputBufferChunkSizes[j] = inputChunkSizes[j + i * num_threads];
                    inputBufferChunkOffsets[j] = inputChunkOffsets[j + i * num_threads];
                }
                else {
                    inputBufferChunkSizes[j] = 0;
                    outputBufferChunkSizes[j] = 0;
                    inputBufferChunkOffsets[j] = 0;
                }
                if (inputBufferChunkSizes[j] > 0) {
                    fseek(inputFile, inputBufferChunkOffsets[j], SEEK_SET);
                    readSize = fread(inputChunk + j * MAX_ENCODED_BUFFER_SIZE, sizeof(unsigned char), inputBufferChunkSizes[j], inputFile);
                    if (readSize != inputBufferChunkSizes[j]) {
                        isDecodingSuccessful = false;
                    }
                }
                omp_unset_lock(&processlock[j]);
            }
            // single last thread as writer
        }
        else if (thread_ID == num_threads + 1) {
            #pragma omp barrier
            #pragma omp single
            // writer
            {
            for (int j = 0;j < num_threads;j++) {
                omp_set_lock(&writelock[j]);
                if (outputBufferChunkSizes[j] > 0) {
                    fwrite(outputChunk + j * (MAX_DECODED_BUFFER_SIZE), sizeof(unsigned char), decodedOutputBufferChunkSizes[j], outputFile);
                }
                omp_unset_lock(&readlock[j]);
            }
            }
        }
        // all other num_threads on processing
        else {
            int work_ID = omp_get_thread_num() - 1;
            #pragma omp critical
            {
                work_ID = current_chunk;
                current_chunk = (current_chunk + 1) % num_threads;
            }
            omp_set_lock(&processlock[work_ID]);
            if (outputBufferChunkSizes[work_ID] > 0) {
                chunkDecoder(inputChunk + work_ID * MAX_ENCODED_BUFFER_SIZE, outputChunk + work_ID * (MAX_DECODED_BUFFER_SIZE), huffmanTree, outputBufferChunkSizes[work_ID], &decodedOutputBufferChunkSizes[work_ID]);
            }
            omp_unset_lock(&writelock[work_ID]);
        }


    }
    for (int j = 0;j < num_threads;j++) {
        omp_destroy_lock(&readlock[j]);
        omp_destroy_lock(&processlock[j]);
        omp_destroy_lock(&writelock[j]);
    }
    free(inputChunk);
    free(outputChunk);
    free(inputBufferChunkSizes);
    free(outputBufferChunkSizes);
    free(decodedOutputBufferChunkSizes);
    free(inputBufferChunkOffsets);
    free(readlock);
    free(processlock);
    free(writelock);
    return isDecodingSuccessful;
}

bool fileDecoderFull(char* inputFileName, char* outputFileName, int num_threads, int mode, int rank) {
    // get byte frequencies in the input file
    FILE* inputFile = fopen(inputFileName, "r");
    if (!inputFile) {
        perror(inputFileName);
        exit(1);
    }
    Dictionary* dict = createDictionary(MAX_HEAP_SIZE);
    getFrequenciesFromEncodedFile(inputFile, dict);

    // get huffman tree and alphabet
    Node* huffmanTree = getHuffmanTree(dict);

    // get chunk offsets
    unsigned long long numChunks;
    unsigned long long* chunkOffsets = getChunkOffsetsFromEncodedFile(inputFile, &numChunks);
    unsigned long long* inputChunkSizes = getOriginalChunkSizesFromEncodedFile(inputFile, &numChunks);

    // prepare input and output files for encoding
    fseek(inputFile, 0, SEEK_SET);
    FILE* outputFile = fopen(outputFileName, "w");
    if (!outputFile) {
        perror(outputFileName);
        exit(1);
    }
    bool isDecodingSuccessful = false;
    if (mode == 0)
        isDecodingSuccessful = fileDecoderBarrier(inputFile, outputFile, huffmanTree, chunkOffsets, inputChunkSizes, numChunks, num_threads);
    else
        isDecodingSuccessful = fileDecoderLocks(inputFile, outputFile, huffmanTree, chunkOffsets, inputChunkSizes, numChunks, num_threads);

    freeNode(huffmanTree);
    huffmanTree = NULL;
    freeDictionary(dict);
    dict = NULL;
    free(chunkOffsets);
    free(inputChunkSizes);
    fclose(inputFile);
    fclose(outputFile);
    return  isDecodingSuccessful;
}

int fileProcesser(char* inputname, char* outputname, int num_threads, bool (*processing)(char* arg1, char* arg2, int num_threads, int mode, int rank), int mode) {
    printf("Processing file %s and saving as %s\n", inputname, outputname);
    bool completed = (*processing)((char*)inputname, (char*)outputname, num_threads, mode, 0);
    if (!completed) {
        printf("Error on file %s\n", outputname);
        exit(1);
    }
    return 0;
}

int main(int argc, char** argv) {
    char inputname[PATH_MAX];
    char outputname[PATH_MAX];
    void* processingFunction = fileDecoderFull;
    int num_threads = 1;
    num_threads = omp_get_max_threads();
    int mode = 0;       // 0 for barrier, 1 for locks
    strcpy(inputname, argv[optind]);
    strcpy(outputname, argv[optind + 1]);
    fileProcesser(inputname, outputname, num_threads, processingFunction, mode);
    return 0;
}