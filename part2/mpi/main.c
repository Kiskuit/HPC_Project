#include "projet.h"
#include <unistd.h>
#include <sys/types.h>

#define RATIO 10
/* TODO estimates how many time we have pruning with a counter in 
 * original. Adjust ratio according to result, speak about it in report */

/* 2017-02-23 : version 1.0 */
/* Global Variables */
unsigned long long int node_searched = 0;

/* Create the MPI_Datatype corresponding to the meta_t struct */
MPI_Datatype *MPI_meta_creator () {
    /* Datatype for tree struct */
    /* Number of different blocks of types in struct */
    const int nbTypes_tree = 3;
    /* Number of elements in each block in struct */
    const int blockLengths_tree[] = {3*128, 11, 1+MAX_DEPTH};
    /* MPI types of each block in struct */
    const MPI_Datatype types_tree[] = {MPI_CHAR, MPI_INT, MPI_UNSIGNED_LONG};
    /* Blocks offsets in struct */
    const MPI_Aint offsets_tree[] = {offsetof (tree_t, pieces),
                                offsetof (tree_t, side),
                                offsetof (tree_t, hash)};
    MPI_Datatype *MPI_tree;
    if ( (MPI_tree = malloc(sizeof(MPI_Datatype))) == NULL) {
        fprintf(stderr, "malloc error in MPI_tree_creator()\n");
        exit(1);
    }
    MPI_Type_create_struct (nbTypes_tree, blockLengths_tree, offsets_tree, types_tree, MPI_tree);
    MPI_Type_commit (MPI_tree);

    /* Datatype for result struct */
    const int nbTypes_res = 1;
    const int blockLengths_res[] = {3+MAX_DEPTH};
    const MPI_Datatype types_res[] = {MPI_INT};
    const MPI_Aint offsets_res[] = {offsetof (result_t, score)};
    MPI_Datatype *MPI_result;
    if ( (MPI_result = malloc(sizeof(MPI_Datatype))) == NULL) {
        fprintf(stderr, "malloc error in MPI_tree_creator()\n");
        exit(1);
    }
    MPI_Type_create_struct (nbTypes_res, blockLengths_res, offsets_res, types_res, MPI_result);
    MPI_Type_commit (MPI_result);

    /* Datatype for meta struct */
    const int nbTypes = 4;
    const int blockLengths[] = {1,1,1,1};
    const MPI_Datatype types[] = {MPI_INT, MPI_UNSIGNED_LONG,
                                    *MPI_tree, *MPI_result};
    const MPI_Aint offsets[] = {offsetof (meta_t, index),
                                offsetof (meta_t, nodes),
                                offsetof (meta_t, tree),
                                offsetof (meta_t, result)};
    MPI_Datatype *MPI_meta;
    if ( (MPI_meta = malloc(sizeof(MPI_Datatype))) == NULL) {
        fprintf(stderr, "malloc error in MPI_tree_creator()\n");
        exit(1);
    }
    MPI_Type_create_struct (nbTypes, blockLengths, offsets, types, MPI_meta);
    MPI_Type_commit (MPI_meta);

    MPI_Type_free (MPI_tree);
    MPI_Type_free (MPI_result);

    return MPI_meta;
}
    

/* Master's core function : TODO EDIT DOC
 * -First : Perform breadth-first search until there are
 * at least 10x more tasks than slaves.
 * -Second : Distribute works to slave as soon as they are ready.
 * -Third : Recombine all work done by slaves.
*/
/* TODO set pruned flag to FALSE everywhere (by default) */

void alpha (tree_t *T, result_t *result) {
    // Leftmost branch is first stored in the tree
    int beg, sizeTree = T->depth+1, branchSize = T->depth+1;
    recTree_t *masterTree;
    if ( (masterTree=malloc(sizeTree * sizeof(recTree_t))) == NULL) {
        fprintf(stderr, "malloc error in alpha()\n");
        exit(1);
    }
    masterTree[0].tree = T;
    masterTree[0].result = result;
    masterTree[0].pruned = FALSE;

    // Array to store moves
    int *n_moves = calloc (sizeTree, sizeof(int));
    move_t **moves = calloc (sizeTree, sizeof(int*));
    if (!n_moves || !moves) {
        fprintf (stderr, "calloc error in alpha()\n");
        exit(1);
    }
    for (int i=0 ; i<sizeTree ; i++) 
        if ( (moves[i]=calloc(MAX_MOVES,sizeof(int))) == NULL) {
            fprintf(stderr, "calloc error in alpha()\n");
            exit(1);
        }


    /* ------------------------------------------------------------ *
     * ---------- Follows the left branch to the bottom ----------- *
     * ------------------------------------------------------------ */
    for (int i=0 ; i<sizeTree ; i++) {
        node_searched++;
        tree_t *tmpTree = masterTree[i].tree;
        result_t *tmpResult = masterTree[i].result;

        masterTree[i].parentId = i-1;

        /* ----- Equivalent of evaluate ----- */
        tmpResult->score = -MAX_SCORE - 1;
        tmpResult->pv_length = 0;
        if (test_draw_or_victory (tmpTree, tmpResult)) {
            branchSize = sizeTree = i+1; /* TODO check that isnt'it just 'i' */
            break;
        }
        if (i==sizeTree-1) {
            tmpResult->score = (2 * tmpTree->side -1) * heuristic_evaluation (tmpTree);
            break;
        }
        compute_attack_squares (tmpTree);
        n_moves[i] = generate_legal_moves (tmpTree, moves[i]);
        if (n_moves[i]==0) {
            tmpResult->score = check(tmpTree)?-MAX_SCORE:CERTAIN_DRAW;
            branchSize = sizeTree = i+1; /* TODO Same above */
            break;
        }
        sort_moves (tmpTree, n_moves[i], mvoes[i]);
        masterTree[i+1].move = moves[i][0];
        masterTree[i+1].tree = malloc (sizeof(tree_t));
        masterTree[i+1].result = malloc (sizeof(result_t));
        masterTree[i+1].pruned = FALSE;
        if (!masterTree[i+1].tree || ! masterTree[i+1].result) {
            fprintf (stderr, "malloc error in alpha()\n");
            exit(1);
        }
        play_move (tmpTree, moves[i][0], masterTree[i+1].tree);
    } // END for (sizeTree)

    // Analyze leftmost leaf
    // TODO simplify notations!
    if (masterTree[sizeTree-1].result->score > masterTree[sizeTree-2].result->score) {
        masterTree[sizeTree-2].result->score = masterTree[sizeTree-1].result->score;
        masterTree[sizeTree-2].result->best_move = masterTree[sizeTree-1].move;
        masterTree[sizeTree-2].result->pv_length = masterTree[sizeTree-1].result->pv_length+1;
        memcpy (masterTree[sizeTree-2].result->PV+1, masterTree[sizeTree-1].result->PV,
                masterTree[sizeTree-1].result->pv_length * sizeof(int));
        masterTree[size-2].result->PV[0] = masterTree[sizeTree-1].move;
    }
    if (masterTree[siezTree-1].result->score > masterTree[sizeTree-2].tree->beta) {
        n_moves[sizeTree-2] = 1;
        masterTree[sizeTree-2].pruned = TRUE;
    }

    /* ------------------------------------------------------------ *
     * ----- Go up the tree, create task, and distribute work ----- *
     * ------------------------------------------------------------ */
    beg = sizeTree;
    for (int i=sizeTree-2 ; i>0 ; i--) {
        // TODO change to pointer to save the copy
        rectTree_t child = masterTree[i], parent = masterTree[i-1];
        int childScore = -child.result->score;
        // TODO check with and without the if
        if (childScore > parent.result->score) {
            parent.result->score = childScore;
            parent.result->best_move = child.move;
            parent.result->pv_length = child.result->pv_length + 1;
            memcpy (parent.result->PV+1, child.result->PV,
                    child.result->pv_length * sizeof(int));
            parent.result->PV[0] = child.move;
        }

        // Prunning
        if (childR->score >= parentT->beta) {
            n_moves[i-1] = 1; // TODO remove this?
            masterTree[i-1].pruned = TRUE;
        }
        // If it has not been pruned...
        if (!masterTree[i].pruned) {
            int lower = sizeTree;
            sizeTree += n_moves[i] - 1;
            nb_tasks += n_moves[i] - 1;
            if ( (masterTree = realloc(masterTree, sizeTree*sizeof(recTree_t))) == NULL) {
                fprintf(stderr, "realloc error in alpha()\n");
                exit(1);
            }
            // ... Add every child but the first to the taskslist
            for (int j=lower ; j<sizeTree ; j++) {
                masterTree[j].parentId = i;
                masterTree[j].move = moves[j-lower+1];
                masterTree[j].tree = malloc (sizeof(tree_t));
                masterTree[j].result = malloc (sizeof(result_t));
                if (!masterTree[j].tree || !masterTree[j].result) {
                    fprintf(stderr, "malloc error in alpha()\n");
                    exit(1);
                }
                play_move(masterTree[i].tree, masterTree[j].move, masterTree[j].tree);
            }
        }
        // If enough tasks, distribute tasks
        if (nb_tasks >= RATIO*nb_proc) {
            distribute_work (beg, sizeTree, masterTree);
            nb_tasks = 0;
            beg = sizeTree;
        }
    }
    /* -------------------------------------------- *
     * ----- Recombine results of left branch ----- *
     * -------------------------------------------- */
    for (int i=branchSize-1 ; i>0 ; i++) {
        recTree_t *child = &masterTree[i], *parent = &masterTree[i-1];
        int childScore = -child->result->score;
        if (childScore > parent->result->score) {
            parent->result->score = childScore;
            parent->result->best_move = child->move;
            parent->result->pv_length = child->result->pv_length+1;
            memcpy (parent->result->PV+1, child->result->PV,
                    child->result->pv_length*sizeof(int));
            parent->result->PV[0] = child->move;
        }
        parent->tree->alpha = MAX(parent->tree->alpha, childScore);
    }
    
    // Memory
    for (int i=0 ; i<T->depth+1 ; i++)
        free(moves[i]);
    free(moves);
    free(n_moves);
}

void distribute_work (int start, int sizeTree, recTree_t *masterTree) {
    MPI_Datatype *MPI_meta_tye = MPI_meta_creator();
    int nb_proc;
    MPI_Comm_size (MPI_COMM_WORLD, &nb_proc);

    for (int dest=1 ; dest<nb_proc ; dest++) {
        meta_t metaSend;
        metaSend.index = start;
        metaSend.nodes = 0;
        metaSend.tree = *masterTree[start].tree;
        metaSend.result = *masterTree[start].result;
        MPI_Send (&metaSend, 1, *MPI_meta_type, dest, TAG_CONTINUE,
                MPI_COMM_WORLD);
        start++;
    }

    int slaveFinished = 0;
    while (slaveFinished < nb_proc-1) {
        MPI_Status status;
        meta_t metaRecv;
        MPI_Recv (&metaRecv, 1, *MPI_meta_type, MPI_ANY_SOURCE,
                MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        int dest = status.MPI_SOURCE;

        *masteTree[metaRecv.index].tree = metaRecv.tree;
        *masterTree[metaRecv.index].result = metaRecv.result;
        node_searched += metaRecv.nodes;

        recTree_t *child = &masterTree[metaRecv.index];
        recTree_t *parent = &masterTree[child->parentId];
        int childScore = -child->result->score;

        if (childScore > parent->result->score) {
            parent->result->score = childScore;
            parent->result->best_move = child->move;
            parent->result->pv_length = child->pv_length + 1;
            memcpy (parent->result->PV+1, child->result->PV,
                    child->result->pv_length*sizeof(int));
            parent->result->PV[0] = child->move;
        }
        if (childScore >= parent->tree->beta) {
            parent->pruned = TRUE;
        }
        parent->tree->alpha = MAX(parent->tree->alpha, childScore);

        // Search an element whose parent hasnt been pruned
        while (masterTree[masterTree[start].parentId].pruned)
            start++;
        if (start<sizeTree) { // More tasks to send
            metaRecv.index = start;
            metaRecv.nodes = 0;
            metaRecv.tree = *masterTree[start].tree;
            metaRecv.result = *masterTree[start].result;
            // TODO can we use Isend?
            MPI_Send (&metaRecv, 1, *MPI_meta_type, dest,
                    TAG_CONTINUE, MPI_COMM_WORLD);
            start++;
        }
        else { // TODO check if we really need to send smth here
            MPI_Send (NULL, 0, *MPI_meta_type, dest, TAG_STOP, MPI_COMM_WORLD);
            slaveFinished++;
        }
    }
    MPI_Type_free (MPI_meta_type);
}




// Only slaves call evaluate
void evaluate(tree_t * T, result_t *result)
{
    node_searched++;

    move_t moves[MAX_MOVES];
    int n_moves;

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


void decide(int rank, int nb_proc, tree_t * T, result_t *result){
    /* Master executes this */
    if (rank==0) {
        for (int depth = 1;; depth++) {
            T->depth = depth;
            T->height = 0;
            T->alpha_start = T->alpha = -MAX_SCORE - 1;
            T->beta = MAX_SCORE + 1;

            printf("=====================================\n");
            master(T, result);

            printf("depth: %d / score: %.2f / best_move : ", T->depth, 0.01 * result->score);
            print_pv(T, result);

            if (DEFINITIVE(result->score)) {
                /* Signal all slaves that they can stop & exit */
                for (int dest=1 ; dest<nb_proc ; dest++)
                    MPI_Send (NULL, 0, MPI_INT, dest, TAG_OVER, MPI_COMM_WORLD);
                break;
            }
        }
    }
    /* Slaves execute this */
    else {
        int isOver = FALSE;
        while (!isOver)
            isOver = slave_function();
    }
}

int slave_function() {
    MPI_Status status;
    MPI_Datatype *MPI_meta_type = MPI_meta_creator ();
    int tag;
    meta_t meta;
    const int dest = 0;
    /* Loop to receive work, execute it, and send it back */
    do {
        MPI_Probe (dest, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        tag = status.MPI_TAG;

        if(tag == TAG_CONTINUE) { /* More work to do */
            MPI_Recv (&meta, 1, *MPI_meta_type, dest, tag,
                    MPI_COMM_WORLD, NULL);

            /* Main job */
            evaluate (&meta.tree, &meta.result);

            /* Send response */
            meta.nodes = node_searched;
            MPI_Send (&meta, 1, *MPI_meta_type, dest, TAG_ANS, MPI_COMM_WORLD);
            node_searched = 0;
        }
        else {
            /* Recv to unlock slave */
            MPI_Recv (NULL, 0, *MPI_meta_type, dest, tag,
                    MPI_COMM_WORLD, NULL);
            MPI_Type_free (MPI_meta_type);
            return (tag==TAG_OVER);
        }
    } while (tag == TAG_CONTINUE);
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

    if (comSize < 2) {
        fprintf(stderr, "This program requires at least two processes!\n");
        exit(1);
    }

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
    }

    /* Start computation */
    decide(rank, comSize, &root, &result);

    if (rank==0) {
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

    /* MPI Finalization */
    MPI_Finalize();
    return 0;
}

