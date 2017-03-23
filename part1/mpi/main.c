#include "projet.h"
#include <mpi.h>

/* 2017-02-23 : version 1.0 */

unsigned long long int node_searched = 0;


/*
 * IDEA :
 *  - Master call pre_evaluate
 *      Reccursively goes as deep as needed (~10x nbProc)
 *      Sets up trees and results structures for slaves
 *      Distribute work to slaves.
 *  - Slaves wait for work
 *      call to evaluate when receiving trees,
 *      signal master upon completion
 *      wait for more
 *          -> last message with special TAG to exit and call MPI_Finalize()
 */

void pre_evaluate(/* TODO */) {
    // TODO handle node_searched increment in this function
    int nb_tasks, n_moves;
    move_t moves[MAX_MOVES];

    tree_t *taskTrees, *parentTrees;
    result_t *taskResults, *parentResults;

    result->score = -MAX_SCORE -1;
    result->pv_length = 0;
    if (test_draw_or_victory(T, result))
        return;

    if (T->depth == 0) {
        result->score = (2*T->side -1)*heuristic_evaluation(T);
        return;
    }

    nb_tasks = n_moves = generate_legal_moves(T, moves);
    if ( (taskTrees=calloc(nb_tasks, sizeof(tree_t)))==NULL) {
        fprintf(stderr,"calloc error in pre_evaluate()\n");
        exit(1);
    }
    if ( (taskResults=calloc(nb_tasks, sizeof(tree_t))) == NULL) {
        fprintf(stderr,"calloc error in pre_evaluate()\n");
        exit(1);
    }
    parentTrees = T;
    parentResults = result;
    parentMoves = moves;
    *parentID=4;
    while (nb_tasks < 10*nb_proc) {
        // TODO take into account case where not enough tasks
        // Go deeper
        int sum = 0, j = 0;
        for (int i = 0 ; i<nb_tasks ; i++) {
            // TODO change T to parent somehow
            //      and moves to parentMoves
            if (i >= parentID[j]) 
                j++;
            play_move(parentTrees[j], moves[i], &taskTree[i]);

            taskResults[i].score = -MAX_SCORE - 1;
            taskResults[i].pv_length = 0;
            if (test_draw_or_victory(&taskTrees[i], &taskResults[i])){
                // What are we doing here?
                continue;
            }
            if (taskTrees[i].depth == 0) {
                taskResults[i].score = (2*taskTree[i].side - 1)*heuristic_evaluation(&taskTree[i]);
                // Same as above
                continue;
            }
            n_moves = generate_legal_moves(&taskTrees[i], /*TODO*/);
            nb_tasks += n_moves;
        }
    }
    // Do parallel thingy
    /* 
     * Master sets up to send tasks to every slave
     * Send & Receive
     * Upon reception, make adjustment (best score thingy)
     */
}

// Only slaves call evaluate
void evaluate(tree_t * T, result_t *result)
{
    node_searched++;

    move_t moves[MAX_MOVES];
    int n_moves;

    /* MPI vars */
    int rank, comSize;
    MPI_Comm_size(MPI_COMM_WORLD, &comSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    result->score = -MAX_SCORE - 1;
    result->pv_length = 0;

    if (test_draw_or_victory(T, result))
        return;

    if (TRANSPOSITION_TABLE && tt_lookup(T, result))     /* la réponse est-elle déjà connue ? */
        return;

    compute_attack_squares(T);

    /* profondeur max atteinte ? si oui, évaluation heuristique */
    if (T->depth == 0) {
        result->score = (2 * T->side - 1) * heuristic_evaluation(T);
        return;
    }

    n_moves = generate_legal_moves(T, &moves[0]);
    fprintf(stderr,"Height : %d\n", T->height);

    /* absence de coups légaux : pat ou mat */
    if (n_moves == 0) {
        result->score = check(T) ? -MAX_SCORE : CERTAIN_DRAW;
        return;
    }

    if (ALPHA_BETA_PRUNING)
        sort_moves(T, n_moves, moves);

    /* évalue récursivement les positions accessibles à partir d'ici */
    for (int i = 0; i < n_moves; i++) {
        tree_t child;
        result_t child_result;

        play_move(T, moves[i], &child);

        evaluate(&child, &child_result);

        int child_score = -child_result.score;

        if (child_score > result->score) { 
            result->score = child_score;
            result->best_move = moves[i];
            result->pv_length = child_result.pv_length + 1;
            for(int j = 0; j < child_result.pv_length; j++)
                result->PV[j+1] = child_result.PV[j];
            result->PV[0] = moves[i];
        }

        if (ALPHA_BETA_PRUNING && child_score >= T->beta)
            break;    

        T->alpha = MAX(T->alpha, child_score);
    }

    if (TRANSPOSITION_TABLE)
        tt_store(T, result);
}


void decide(tree_t * T, result_t *result){
        for (int depth = 1;; depth++) {
            T->depth = depth;
            T->height = 0;
            T->alpha_start = T->alpha = -MAX_SCORE - 1;
            T->beta = MAX_SCORE + 1;

            printf("=====================================\n");
            evaluate(T, result);

            printf("depth: %d / score: %.2f / best_move : ", T->depth, 0.01 * result->score);
            print_pv(T, result);

            if (DEFINITIVE(result->score))
                break;
        }

}

int main(int argc, char **argv)
{  
    tree_t root;
    result_t result;
    int rank, comSize;
    /* MPI Initialization */
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    if (rank==0) {
        if (argc < 2) {
            printf("usage: %s \"4k//4K/4P w\" (or any position in FEN)\n", argv[0]);
            exit(1);
        }

        if (ALPHA_BETA_PRUNING)
            printf("Alpha-beta pruning ENABLED\n");

        if (TRANSPOSITION_TABLE) {
            printf("Transposition table ENABLED\n");
            init_tt();
        }

        parse_FEN(argv[1], &root);
        print_position(&root);


        decide(&root, &result);

        printf("\nDécision de la position: ");
        switch(result.score * (2*root.side - 1)) {
            case MAX_SCORE: printf("blanc gagne\n"); break;
            case CERTAIN_DRAW: printf("partie nulle\n"); break;
            case -MAX_SCORE: printf("noir gagne\n"); break;
            default: printf("BUG\n");
        }

        printf("Node searched: %llu\n", node_searched);

        if (TRANSPOSITION_TABLE)
            free_tt();
    }
    else {
        // TODO
        fprintf(stderr,"Proc %d ne fait rien\n", rank);
    }

    /* MPI Finalization */
    MPI_Finalize();
    return 0;
}
