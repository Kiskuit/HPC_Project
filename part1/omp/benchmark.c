#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char** argv) {
    if(argc != 4) {
        printf("%d arguments\n", argc);
        printf("Usage : benchmark nbCalls param nbThreads\n");
        return EXIT_FAILURE;
    }
    int nbTimes = atoi(argv[1]);
    int nbThreads = atoi(argv[3]);
    FILE* fp;
    char output[256];
    char command[256];
    char envSave[256]; 
    char envPara[256];
    char envSeq[256];
    float avParallel=0;
    float avSequentiel=0;
    strcpy(envSave, "OMP_NUM_THREADS=");
    strcat(envSave, getenv("OMP_NUM_THREADS"));
    strcpy(envPara, "OMP_NUM_THREADS=");
    strcat(envPara, argv[3]);
    strcpy(envSeq, "OMP_NUM_THREADS=1");
    strcpy(command, "/usr/bin/time --format='%e' ./decide \"");
    strcat(command, argv[2]);
    strcat(command, "\" 2>&1 > /dev/null");
    printf("command :'%s'\n", command);
    printf("Testing %d times\n\n", nbTimes);
    
    //debut par test
    putenv(envPara);
    for (int i=0; i<nbTimes; i++) {
        fp = popen("/usr/bin/time --format='%e' ./decide \"4k//4K/4P w\" 2>&1 > /dev/null", "r");
        if (fp == NULL) {
            printf("Failed to run command\n" );
            exit(1);
        }

        if(fgets(output, sizeof(output)-1, fp) == NULL) return EXIT_FAILURE;
        float value = atof(output);
        avParallel+=value;
    }
    avParallel/=nbTimes;
    
    printf("Average parallel time = %f\n\n", avParallel);
    
    //debut seq test
    putenv(envSeq);
    for (int i=0; i<nbTimes; i++) {
        fp = popen("/usr/bin/time --format='%e' ./decide \"4k//4K/4P w\" 2>&1 > /dev/null", "r");
        if (fp == NULL) {
            printf("Failed to run command\n" );
            exit(1);
        }

        if(fgets(output, sizeof(output)-1, fp) == NULL) return EXIT_FAILURE;
        float value = atof(output);
        avSequentiel+=value;
    }
    avSequentiel/=nbTimes;

    printf("Average seq time = %f\n\n", avSequentiel);
    
    printf("Taux = %f\n", (avSequentiel / avParallel) / nbThreads);
    pclose(fp);
    putenv(envSave);

    return EXIT_SUCCESS;
}
