#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <malloc.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <pthread.h>
#include <errno.h>

struct Info {
    int N, i, chankSize;
    pthread_mutex_t * mutexCounter;
    int *counter;
    double* arr1;
    double* arr2;
    int* start;
};

struct Compare { 
    double val; 
    int index; 
};

double get_wtime() {
    struct timeval T;
    gettimeofday(&T, NULL);
    return 1000 * T.tv_sec + T.tv_usec / 1000;
}

void printArr(double* arr, int size, char* str) {
    int i = 0;
    printf("%s:", str);
    for(i = 0; i < size; i ++) {
        printf(" %f ", arr[i]);
    }
    printf("\n");
}

void printArr2(int* arr, int size, char* str) {
    int i = 0;
    printf("%s:", str);
    for(i = 0; i < size; i ++) {
        printf(" %d ", arr[i]);
    }
    printf("\n");
}

void * printRunComplite(void * CountItem) {
    long TAll = get_wtime(); 
    printf("Program complite on %f%%\n", ((int*)CountItem)[0] * 100.0 / ((int*)CountItem)[1]);
    while(((int*)CountItem)[0] < ((int*)CountItem)[1]) {
        if(get_wtime() - TAll > 1000) {
            printf("Program complite on %f%%\n", ((int*)CountItem)[0] * 100.0 / ((int*)CountItem)[1]);
            TAll = get_wtime();
        }
    }
    return 0;
}

void Sort(int start, int end, double* M2, int* CountItem, pthread_mutex_t * mutexCounter) {
    struct Compare max;
    int j, k;
    double temp = 0;
    pthread_mutex_lock(mutexCounter);
    CountItem[0]++;
    pthread_mutex_unlock(mutexCounter);
    for (j = start; j > end; --j) {
        max.val = M2[j];
        max.index = j;
        for (k = j - 1; k >= end; --k) {
            if (M2[k] > max.val) {
                max.val = M2[k];
                max.index = k;
            }
        }
        temp = M2[j];
        M2[j] = max.val;
        M2[max.index] = temp;
        pthread_mutex_lock(mutexCounter);
        CountItem[0]++;
        pthread_mutex_unlock(mutexCounter);
    }
}
// Этапы
void* generateM1(void * arg) {
    struct Info info = *(struct Info*)arg;
    int i = 0, j = 0, buffStart = 0, buffStop= 0;
    unsigned int randNum = 0;
    while(info.start[i] != -1) {
        buffStart = info.start[i];
        buffStop = (info.start[i] + info.chankSize > info.N)? info.N: info.start[i] + info.chankSize;

        for(j = buffStart; j < buffStop; j ++) {
            randNum = j + info.i * (info.N * 1.5);
            info.arr1[j] = (rand_r(&randNum) % 27000)/100.0 + 1;
            pthread_mutex_lock(info.mutexCounter);
            *info.counter = *info.counter + 1;
            pthread_mutex_unlock(info.mutexCounter);
        }
        
        i++;
    }
    pthread_exit(NULL);
}

void* generateM2(void * arg) {
    struct Info info = *(struct Info*)arg;
    int i = 0, j = 0, buffStart = 0, buffStop= 0;
    unsigned int randNum = 0;
    while(info.start[i] != -1) {
        buffStart = info.start[i] + info.N;
        buffStop = info.N + ((info.start[i] + info.chankSize > info.N / 2)? info.N / 2: info.start[i] + info.chankSize);
        
        for(j = buffStart; j < buffStop; j ++) {
            randNum = j + info.i * (info.N * 1.5);
            info.arr1[j - info.N] = (rand_r(&randNum) % 243000)/100.0 + 270;
            info.arr2[j - info.N] = info.arr1[j - info.N];
            pthread_mutex_lock(info.mutexCounter);
            *info.counter = *info.counter + 1;
            pthread_mutex_unlock(info.mutexCounter);
        }
        
        i++;
    }
    pthread_exit(NULL);
}

void * MapM1(void * arg) {
    struct Info info = *(struct Info*)arg;
    int i = 0, j = 0, buffStart = 0, buffStop= 0;
    //printArr2(info.start, 4, "Start Arr");
    while(info.start[i] != -1) {
        buffStart = info.start[i];
        buffStop = (info.start[i] + info.chankSize > info.N)? info.N: info.start[i] + info.chankSize;
        for(j = buffStart; j < buffStop; j++) {
            info.arr1[j] = sqrt(sinh(info.arr1[j]));
            pthread_mutex_lock(info.mutexCounter);
            *info.counter = *info.counter + 1;
            pthread_mutex_unlock(info.mutexCounter);
        }
        i++;
    }
    pthread_exit(NULL);
}

void * MapM2(void * arg) {
    struct Info info = *(struct Info*)arg;
    int i = 0, j = 0, buffStart = 0, buffStop= 0;
    while(info.start[i] != -1) {
        buffStart = info.start[i];
        buffStop = (info.start[i] + info.chankSize > info.N/2)? info.N/2: info.start[i] + info.chankSize;
        if(buffStart == 0){
            buffStart++;
        }
        for(j = buffStart; j < buffStop; j++) {
            info.arr1[j] = fabs(cos(info.arr2[j] + info.arr2[j - 1]) / sin(info.arr2[j] + info.arr2[j - 1]));
            pthread_mutex_lock(info.mutexCounter);
            *info.counter = *info.counter + 1;
            pthread_mutex_unlock(info.mutexCounter);
        }
        i++;
    }
    pthread_exit(NULL);
}

void * Merge(void * arg) {
    struct Info info = *(struct Info*)arg;
    int i = 0, j = 0, buffStart = 0, buffStop= 0;

    while(info.start[i] != -1) {
        buffStart = info.start[i];
        buffStop = (info.start[i] + info.chankSize > info.N/2)? info.N/2: info.start[i] + info.chankSize;
        if(buffStart == 0){
            buffStart++;
        }
        for(j = buffStart; j < buffStop; j++) {
            info.arr2[j] = fabs(info.arr1[j] - info.arr2[j]);
            pthread_mutex_lock(info.mutexCounter);
            *info.counter = *info.counter + 1;
            pthread_mutex_unlock(info.mutexCounter);
        }
        i++;
    }
    pthread_exit(NULL);
}

void * mainSort(void * arg) {
    struct Info info = *(struct Info*)arg;
    int i = 0, j = 0, buffStart = 0, buffStop= 0;
    int startArr[info.i], endArr[info.i];

    while(info.start[i] != -1) {
        buffStart = info.start[i];
        buffStop = (info.start[i] + info.chankSize > info.i)? info.i: info.start[i] + info.chankSize;
        for(j = buffStart; j < buffStop; j++) {
            if(j == info.i - 1) {
                startArr[j] = info.N/2 - (info.N/2 / info.i * j) - 1;
                endArr[j] = 0;
            }
            else {
                startArr[j] = info.N/2 - (info.N/2 / info.i * j) - 1;
                endArr[j] = info.N/2 - (info.N/2 / info.i * (j + 1));
            }
            Sort(startArr[j], endArr[j], info.arr1, info.counter, info.mutexCounter);
        }
        i++;
    }

    pthread_exit(NULL);
}

void * Reduce(void * arg) {
    struct Info info = *(struct Info*)arg;
    int i = 0, j = 0, buffStart = 0, buffStop= 0;
    double * X = (double*)malloc(sizeof(double) * 1);
    *X = 0; 
    while(info.start[i] != -1) {
        buffStart = info.start[i];
        buffStop = (info.start[i] + info.chankSize > info.N/2)? info.N/2: info.start[i] + info.chankSize;
        for(j = buffStart; j < buffStop; j++) {
            if(((int)(info.arr1[j] / info.arr2[0])) % 2 == 0) {
                *X += sin(info.arr1[j]);
            }
            pthread_mutex_lock(info.mutexCounter);
            *info.counter = *info.counter + 1;
            pthread_mutex_unlock(info.mutexCounter);
        }
        i++;
    }
    pthread_exit((void*)X);
}

// Очередь
int* staticOrder(int chankSize, int maxThread, int curThread, int curPos, int maxPos) {
    int i = 0;
    int* begins;

    if(chankSize == 0) {
        begins = (int*)malloc(sizeof(int) * 1);
        begins[0] = maxPos / maxThread;
        return begins;
    }

    begins = (int*)malloc((maxPos / chankSize + 2) * sizeof(int));
    if(curPos != 0 || curThread !=0) {
        if(curPos == -1 || curPos >= chankSize * curThread) {
            begins[0] = -1;
            return begins;
        }
    }
    curPos = 0;
    while(1) {
        curPos += curThread * chankSize;
        if(curPos >= maxPos) {
            break;
        }
        begins[i] = curPos;
        i++;
        curPos += (maxThread - curThread) * chankSize;
    }
    begins[i] = -1;
    return begins;
}

int* dinamicOrder(int chankSize, int maxThread, int curThread, int curPos, int maxPos) {
    int* begins = (int*)malloc(2 * sizeof(int));

    if(chankSize == 0) {
        begins[0] = 1;
    }
    else {
        if(curPos == -1 || curPos >= maxPos) {
            begins[0] = -1;
        }
        else if(curThread == 0 && curPos == 0) {
            begins[0] = 0;
        }
        else if(curPos + chankSize < maxPos) {
            begins[0] = curPos + chankSize;
        }
    }
    begins[1] = -1;

    return begins;
}

double LoadProcess(void* (*action)(void*), int* (*order)(int, int, int, int, int), int lenght, int countThread, int chankSize, struct Info baseInfo, int isReturn) {
    int curPos = 0, i = 0, prov = 1; //buffStatus;
    double res = 0;
    double* X = 0;
    struct Info info[countThread];
    pthread_t threads[countThread];

    if(chankSize == 0) {
        int * buff = order(0, countThread, 0, 0, lenght);
        chankSize = buff[0];
        free(buff);
    }

    for(i = 0; i < countThread; i ++) {
        info[i].arr1 = baseInfo.arr1;
        info[i].arr2 = baseInfo.arr2;
        info[i].N = baseInfo.N;
        info[i].i = baseInfo.i;
        info[i].counter = baseInfo.counter;
        info[i].mutexCounter = baseInfo.mutexCounter;
        info[i].chankSize = chankSize;
        threads[i] = -1;
    }
    
    prov = countThread;
    while(prov != 0) {
        for(i = 0; i < countThread; i ++) {
            if(threads[i] == -1) {
                info[i].start = order(chankSize, countThread, i, curPos, lenght);
                curPos = info[i].start[0];
                pthread_create(&threads[i], NULL, action, &info[i]);
            }
            else if(pthread_tryjoin_np(threads[i], (void*)&X) != EBUSY) {
                prov --;
                free(info[i].start);

                if(isReturn) {
                    res += *X;
                    free(X);
                }

                info[i].start = order(chankSize, countThread, i, curPos, lenght);
                curPos = info[i].start[0];

                if(info[i].start[0] != -1) {
                    prov ++;
                    pthread_create(&threads[i], NULL, action, &info[i]);
                }
            }
        }
    }
    return res;
}

struct Info createInfo(double arr1[], double arr2[], pthread_mutex_t * mutexCounter, int N, int i, int *counter) {
    struct Info baseInfo;
    baseInfo.arr1 = arr1;
    baseInfo.arr2 = arr2;
    baseInfo.N = N;
    baseInfo.i = i;
    baseInfo.counter = counter;
    baseInfo.mutexCounter = mutexCounter;
    return baseInfo;
}

int main(int argc, char* argv[])  {
    int i = 50, j, k, N, countCore = 1, curIndex = 0, chankSize = 2;
    double X = 0, temp = 0, T1, T2, T3, T4;
    long delta_ms = 0, generateTime = 0, mapTime = 0, mergeTime = 0, sortTime = 0, reduceTime = 0;
    pthread_t threadCounter;
    pthread_mutex_t mutexCounter;
    int CountItem[2];

    struct Info baseInfo;

    N = atoi(argv[1]);
    countCore = atoi(argv[2]);
    chankSize = atoi(argv[3]);

    CountItem[0] = 0;
    CountItem[1] = 5 * N * i;
    
    double M1[N], M2[N/2], M3[N/2];
    int startArr[countCore], endArr[countCore];
    
    pthread_mutex_init(&mutexCounter, NULL);
    pthread_create(&threadCounter, NULL, printRunComplite, CountItem);

    T1 = get_wtime();

    for (i = 0; i < 50; i ++) {
        T3 = get_wtime();
        //Generate
        baseInfo = createInfo(M1, NULL, &mutexCounter, N, i, &CountItem[0]);
        LoadProcess(generateM1, dinamicOrder, N, countCore, 100, baseInfo, 0);
        baseInfo = createInfo(M2, M3, &mutexCounter, N, i, &CountItem[0]);
        LoadProcess(generateM2, dinamicOrder, N/2, countCore, 100, baseInfo, 0);

        T4 = get_wtime();
        generateTime += T4 - T3;
/*
        printf("%d.Generate:\n", i);
        printArr(M1, N, "M1: ");
        printArr(M2, N/2, "M2: ");
        printArr(M3, N/2, "M3: ");
*/
        //Map 1 + ((270 mod 47) mod 7) = 1; 1 + ((270 mod 47) mod 8) = 4
        T3 = get_wtime();
        baseInfo = createInfo(M1, NULL, &mutexCounter, N, i, &CountItem[0]);
        LoadProcess(MapM1, staticOrder, N, countCore, chankSize, baseInfo, 0);

        M2[0] = fabs(cos(M3[0]) / sin(M3[0]));
        CountItem[0]++;

        baseInfo = createInfo(M2, M3, &mutexCounter, N, i, &CountItem[0]);
        LoadProcess(MapM2, staticOrder, N/2, countCore, chankSize, baseInfo, 0);

        T4 = get_wtime();
        mapTime += T4 - T3;
/*
        printf("%d.Map:\n", i);
        printArr(M1, N, "M1: ");
        printArr(M2, N/2, "M2: ");
*/
        //Merge 1 + ((270 mod 47) mod 6) = 6;
        T3 = get_wtime();
        CountItem[0]++;
        baseInfo = createInfo(M1, M2, &mutexCounter, N, i, &CountItem[0]);
        LoadProcess(Merge, staticOrder, N/2, countCore, chankSize, baseInfo, 0);

        T4 = get_wtime();
        mergeTime += T4 - T3;
/*
        printf("%d.Merge:\n", i);
        printArr(M2, N/2, "M2: ");
*/
        //Sort
        T3 = get_wtime();
        for(j = 0; j < countCore; j ++) {
            if(j == countCore - 1) {
                startArr[j] = N/2 - (N/2 / countCore * j) - 1;
                endArr[j] = 0;
            }
            else {
                startArr[j] = N/2 - (N/2 / countCore * j) - 1;
                endArr[j] = N/2 - (N/2 / countCore * (j + 1));
            }
        }

        baseInfo = createInfo(M2, NULL, &mutexCounter, N, countCore, &CountItem[0]);
        LoadProcess(mainSort, staticOrder, countCore, countCore, 1, baseInfo, 0);

        for(j = 0; j < N/2; j ++) {
            temp = DBL_MAX;
            for(k = 0; k < countCore; k ++){
                if(endArr[k] <= startArr[k] && M2[endArr[k]] < temp) {
                    curIndex = k;
                    temp = M2[endArr[k]];
                }
            }
            endArr[curIndex] ++;
            M3[j] = temp;
            CountItem[0]++;
        }
        T4 = get_wtime();
        sortTime += T4 - T3;
/*
        printf("%d.Sort:\n", i);
        printArr(M2, N/2, "M2: ");
        printArr(M3, N/2, "M3: ");
*/
        //Reduce
        T3 = get_wtime();
        for (j = 0; j < N/2 - 1; j++)  {
            if(M3[j] != 0) {
                temp = M3[j];
                break;
            }
            j++;
        }

        baseInfo = createInfo(M3, &temp, &mutexCounter, N, i, &CountItem[0]);
        X += LoadProcess(Reduce, staticOrder, N/2, countCore, chankSize, baseInfo, 1);

        T4 = get_wtime();
        reduceTime += T4 - T3;
/*
        printf("%d.Reduce: %f\n", i + 1, X);
*/
    }
    
    T2 = get_wtime();
    delta_ms =  T2 - T1;
    printf("\nN=%d. All: %ld Gen: %ld Map: %ld Merge: %ld Sort: %ld Red: %ld X=%f\n", N, delta_ms, generateTime, mapTime, mergeTime, sortTime, reduceTime, X);
    pthread_mutex_destroy(&mutexCounter);
    pthread_exit(NULL);
}
