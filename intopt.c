#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define EPSILON 1e-6

//348 lines executed when running GCC -fprofile-arcs (97.7% coverage)

typedef struct simplex_t 
{
    int         m;      // constraints
    int         n;      // decision variables
    int*        var;    // [n+m+1] nonbasic
    double**    a;      // [m][n+1] matrix
    double*     b;      // [m] b array
    double*     x;      // [n+1] x vector
    double*     c;      // [n] coefficients
    double      y;      // y
} simplex_t;

typedef struct node_t 
{
    int         m;      // constraints
    int         n;      // decision variables
    int         k;      // Parent branches on xk.
    int         h;      // branch on xh
    double      xh;     // xh
    double      ak;     // parent ak
    double      bk;     // parent bk
    double*     min;    // [n] lower bounds
    double*     max;    // [n] upper bounds
    double**    a;      // [m][n+1] matrix
    double*     b;      // [m] b array
    double*     x;      // [n] x vector
    double*     c;      // [n] coefficients
    double      z;      // z solution
} node_t;

typedef struct Node 
{
    struct node_t *data;   // pointer to node data
    struct Node *next;  // pointer to next node in list
} Node;


/**
 * @brief removes and returns the first node's data from the linked list
 * 
 * the pop function removes the head of the linked list, obtains the data
 * and updates the head pointer so that it points to the next node in the
 * list. To avoid memory leaks the node that has been removed is freed.
 * 
 * @param head pointer to the pointer of the head of the linked list
 * @return pointer to the data stored in the popped node.
 */
node_t* pop(Node** head) 
{
    Node *removed         = *head;
    *head                 = removed->next;
    node_t *removedData   = removed->data;
    
    free(removed);
    return removedData;
    
}

/**
 * @brief pushes a new node at the beginning of a linked list
 * 
 * push function creates a new node with the given data parameter and inserts
 * it at the beginning of the linked list. If the list is empty, the new node
 * becomes the head of the linked list.
 * 
 * @param head pointer to the pointer of the head of the linked list.
 * @param data pointer to the data to be stored in the new node
 */
void push(Node** head, node_t* data) 
{
    Node *newNode    = calloc(1, sizeof(Node));

    newNode->data    = data;

    if (*head != NULL)
    {
        newNode->next = *head;
    }

    *head = newNode;
}

/* ============================================================================================================================================================================
 *                                                                          Function Prototypes
 * ============================================================================================================================================================================
 */
/* Simplex functions*/
double  simplex(int m, int n, double** a, double* b, double* c, double* x, double y);
double  xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h);
int     init(simplex_t* s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var);
int     pivot(simplex_t* s, int row, int col);
int     initial(simplex_t* s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var);
int     select_nonbasic(simplex_t* s);
void    prepare(simplex_t* s, int k);

/* Branch and bound functions*/
node_t* extend(node_t* p, int m, int n, double** a, double* b, double* c, int k, double ak, double bk);
node_t* initial_node(int m, int n, double** a, double* b, double* c);
void    bound(node_t* p, Node** h, double* zp, double* x);
void    succ (node_t* p, Node** h, int m, int n, double** a, double* b, double* c, int k, double ak, double bk, double* zp, double* x);
double  intopt(int m, int n, double** a, double* b, double* c, double* x);
int     is_integer(double* xp);
int     integer(node_t* p);

/* ============================================================================================================================================================================
 *              
 * ============================================================================================================================================================================
 */

/**
 * @brief initalizes the simplex struct with the given parameters
 * 
 * @param s Pointer to the simplex structure to initialize.
 * 
 * @see simplex_t struct
 * 
 * @return The index `k` of the constraint with the smallest b value
 */
int init(simplex_t* s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var)
{
    int i;
    int k;

    s->m = m;
    s->n = n;
    s->a = a;
    s->b = b;
    s->c = c;
    s->x = x;
    s->y = y;
    s->var = var;

    if(s->var == NULL)
    {
        s->var = calloc(m + n + 1, sizeof(int));

        for (i = 0; i < m + n; i++)
        {
            s->var[i] = i;
        }
    }

    for (k = 0, i = 1; i < m; i++)
    {
        if (s->b[i] < s->b[k])
        {
            k = i;
        }
    }
    return k;
}

/**
 * @brief initializes the simplex tableau and checks solution feasibility
 * 
 * prepares the simplex tableau to solve the linear programming problem.
 * It does this by handling aux variables, checking feasiblity, and updating
 * the tablue to remove variabels if needed
 * 
 * @param s pointer to simplex structure
 * @see simplex_t struct
 * 
 * @return 1 if the problem is feasible, otherwise 0.
 */
int initial(simplex_t* s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var)
{
    int i    = 0;
    int j    = 0;
    int k    = 0;
    double w = 0.0;

    k = init(s, m, n, a, b, c, x, y, var);

    if(b[k] >= 0){ return 1; }

    prepare(s, k);

    n = s->n;

    s->y = xsimplex(m, n, s->a, s->b, s->c, s->x, 0.0, s->var, 1);

    for(i = 0; i < m + n; i++)
    {
        if(s->var[i] == m + n - 1)
        {
            if(fabs(s->x[i]) > EPSILON)
            {
                free(s->x);
                free(s->c);
                s->x = NULL;
                s->c = NULL;
                return 0;
            } else
            {
                break;
            }
        }
    }

    if (i >= n)
    {
        for (j = k = 0; k < n; k++)
        {
            if(fabs(s->a[i - n][k]) > fabs(s->a[i - n][j]))
            {
                j = k;
            }
            
        }

        pivot(s, i - n, j);
        i = j;
    }

    if(i < n-1)
    {
        k                   = s->var[i];
        s->var[i]           = s->var[n - 1]; 
        s->var[n - 1]       = k;

        for(k = 0; k < m; k++)
        {
            w               = s->a[k][n - 1];
            s->a[k][n - 1]  = s->a[k][i];
            s->a[k][i]      = w;
        }

    }

    free(s->c);
    s->c = c;
    s->y = y;

    for (k = n - 1; k < n + m - 1; k++)
    {
        s->var[k] = s->var[k+1];
    }

    n = s->n = s->n - 1;
    double* t = calloc (n, sizeof(double));

    for (k = 0; k < n; k++)
    {
        for (j = 0; j < n; j++)
        {
            if(k == s->var[j])
            {
                t[j] = t[j] + s->c[k];
                goto next_k;
            }
        }

        for (j = 0; j < m; j++)
        {
            if(s->var[n + j] == k)
            {
                break;
            }
        }

        s->y = s->y + s->c[k] * s->b[j];

        for (i = 0; i < n; i++)
        {
            t[i] = t[i] - s->c[k] * s->a[j][i];
        }

        next_k:;

    }

    for (i = 0; i < n; i++)
    {
        s->c[i] = t[i];
    }

    free(t);
    free(s->x);
    t = NULL;
    s->x = NULL;

    return 1;
}

/**
 * @brief selects a non-basic variable.
 * 
 * iterates through coefficients to find first positive non-basic variable
 * 
 * @param s pointer to the simplex struct
 * @return index of the first positive non-basic variable
 */
int select_nonbasic(simplex_t* s)
{
    for(int i = 0; i < s->n; i++)
    {
        if(s->c[i] > EPSILON)
        {
            return i;
        }
    }
    return -1;
}

/**
 * @brief perfoms pivot operation for the simplex tableau
 * 
 * updates the simplex tableau after selecting pivot row and column
 * 
 * @param s pointer to the simplex struct
 * @param row index of the pivot row
 * @param col index of the pivot col
 */
int pivot(simplex_t* s, int row, int col)
{
    double** a  = s->a;
    double* b   = s->b;
    double* c   = s->c;
    int m       = s->m;
    int n       = s->n;
    int i       = 0;
    int j       = 0;
    int t       = 0;

    t               = s->var[col];
    s->var[col]     = s->var[n+row];
    s->var[n+row]   = t;

    s->y = s->y + c[col] * b[row] / a[row][col];

    for (i = 0; i < n; i++)
    {
        if(i != col)
        {
            c[i] -= c[col] * a[row][i] / a[row][col];
        }
    }

    c[col] = -c[col] / a[row][col];

    for (i = 0; i < m; i++)
    {
        if(i != row)
        {
            b[i] -= a[i][col] * b[row] / a[row][col];
        }
    }

    for(i = 0; i < m; i++)
    {
        if(i != row)
        {
            for (j = 0; j < n; j++)
            {
                if(j != col)
                {
                    a[i][j] -= a[i][col] * a[row][j] / a[row][col];
                }
            }
        }
    }

    for (i = 0; i < m; i++)
    {
        if(i != row)
        {
            a[i][col] = -a[i][col] / a[row][col];
        }
    }

    for (i = 0; i < n; i++)
    {
        if(i != col)
        {
            a[row][i] = a[row][i] / a[row][col];

        }
    }
    
    b[row] = b[row] / a[row][col];
    a[row][col] = 1 / a[row][col];
}

/**
 * @brief prepares the simplex tableau for handling artifical variables
 * 
 * adds aux. variable to the simplex tableau in order to ensure feasiblity
 * of solution.
 * 
 * @param s pointer to the simplex struct
 * @param k index of the constraint with smallest b value.
 */
void prepare(simplex_t* s, int k)
{
    int m = s->m;
    int n = s->n;

    int i = 0;

    for (i = m + n; i > n; i--)
    {
        s->var[i] = s->var[i - 1];
    }

    s->var[n] = m + n;
    n++;

    for (int i = 0; i < m; i++)
    {
        s->a[i][n-1 ] = -1.0;
    }

    s->x = calloc(m + n, sizeof(double));
    s->c = calloc(n,     sizeof(double));

    s->c[n-1]   = -1.0;
    s->n        = n;

    pivot(s, k, n - 1);
}

/**
 * @brief solves the linear program problem using the simplex method
 * 
 * xsimples handles the aux. variables and feasibility through pivot
 * operations iteratively  until the optimal value is obtained.
 * Either the optimal solution is found, or its unbounded or its infeasible
 * 
 * @param h
 * @see simplex struct
 * 
 * @return Optimal value of the objective function, INFINITY if unbounded, or NaN if infeasible
 */
double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h)
{
    simplex_t s;
    int i, row, col;

    if (!initial(&s, m, n, a, b,c, x, y, var))
    {
        free(s.var);
        s.var = NULL;
        return NAN;
    }

    while ((col = select_nonbasic(&s)) >= 0)
    {
        row = -1;

        for (i = 0; i < m; i++)
        {
            if (a[i][col] > EPSILON && (row < 0 || b[i] / a[i][col] < b[row] / a[row][col]))
            {
                row = i;
            }
        }

        if (row < 0)
        {
            free(s.var);
            s.var = NULL;
            return INFINITY;
        }

        pivot (&s, row, col);
    }

    if (h == 0){
        for (i = 0; i < n; i++)
        {
            if (s.var[i] < n)
            {
                x[s.var[i]] = 0;
            }
        }

        for (i = 0; i < m; i++)
        {
            if (s.var[n + i] < n)
            {
                x[s.var[n + i]] = s.b[i];
            }   
        }
        free(s.var);
        s.var = NULL;
                         
    } else
    {
        for (i = 0; i < n; i++)
        {
            x[i] = 0;
        }

        for (i = n; i < n + m; i++)
        {
            x[i] = s.b[i - n];
        }
    }
        
    return s.y;
}

/**
 * @brief solves the linear programming problem using the simplex method
 * 
 * this function is a wrapper function for the xsimplex function
 * 
 * @see xsimplex
 * 
 * @return optimal value of the objective function
 */
double simplex(int m, int n, double** a, double* b, double* c, double* x, double y)
{
    return xsimplex(m, n, a , b, c, x, y, NULL, 0);
}

node_t* initial_node(int m, int n, double** a, double* b, double* c){

    node_t* p   = calloc(1, sizeof(node_t));
    p->a        = calloc(m + 1, sizeof(double *));

    for (int i = 0; i < m + 1; i++)
    {
        p->a[i] = calloc(n + 1, sizeof(double));
    }
    
    p->b    = calloc(m+1,   sizeof(double));
    p->c    = calloc(n+1,   sizeof(double));
    p->x    = calloc(n+1,   sizeof(double));
    p->min  = calloc(n,     sizeof(double));
    p->max  = calloc(n,     sizeof(double));

    p->m    = m;
    p->n    = n;

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            p->a[i][j] = a[i][j];
        }
        
    }
    
    memcpy(p->b, b, sizeof(double) * m);
    memcpy(p->c, c, sizeof(double) * n);
    

    for(int i = 0; i < n; i++)
    {
        p->min[i] = -INFINITY;
        p->max[i] = INFINITY;
    }

    return p;
}

/**
 * @brief extends an node with more constraints and bounds
 * 
 * the function creates and allocates memory for a new node, copying the 
 * constraints, bounds and objective coefficient from node p. It then adds
 * a new constraint defined by the parameters ak and bk. The function also
 * updates the max and min bounds and adjusts the constraint matrix b values.
 * 
 * @param p pointer to parent node
 * @see simplex struct and node struct
 * 
 * @return pointer to the extended node.
 */
node_t* extend(node_t* p, int m, int n, double** a, double* b, double* c, int k, double ak, double bk)
{
    node_t* q = calloc(1, sizeof(node_t));
    int i,j;

    q->k      = k;
    q->ak     = ak;
    q->bk     = bk;

    if (ak > 0 && p->max[k] < INFINITY) { q->m = p->m; }
    else if (ak < 0 && p->min[k] > 0)   { q->m = p->m; }
    else                                { q->m = p->m + 1; }

    q->n = p->n;
    q->h = -1;

    q->a    = calloc(q->m + 1,  sizeof(double*));
    q->b    = calloc(q->m + 1,  sizeof(double));
    q->c    = calloc(q->n + 1,  sizeof(double));
    q->x    = calloc(q->n + 1,  sizeof(double));
    q->min  = calloc(n,         sizeof(double));
    q->max  = calloc(n,         sizeof(double));

    for (int i = 0; i < q->m + 1; i++)
    {
        q->a[i] = calloc(q->n + 1, sizeof(double));
    }

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            q->a[i][j] = a[i][j];
        }
        q->b[i] = b[i];
    }

    for (int i = 0; i < n; i++)
    {
        q->c[i] = c[i];
    }
    
    for (int i = 0; i < p->n; i++)
    {
        q->min[i] = p->min[i];
        q->max[i] = p->max[i];
    }

    
    if (ak > 0)
    {
        if (q->max[k] == INFINITY || bk < q->max[k])
        {
            q->max[k] = bk;
        }
    } else if (q->min[k] == -INFINITY || -bk > q->min[k])
    {
        q->min[k] = -bk;
    }

    for (i = m, j = 0; j < n; j++)
    {
        if (q->min[j] > -INFINITY)
        {
            q->a[i][j] = -1;
            q->b[i] = -q->min[j];
            i++;
        }
        if (q->max[j] < INFINITY) {
            q->a[i][j] = 1;
            q->b[i] = q->max[j];
            i++;
        }
    }
    return q;
}

/**
 * @brief checks if a value is an integer by comparing the input value with
 * its nearest round integer. Uses user defined epsilon 1e-6
 * 
 * @param xp pointer to the double value that is to be checked.
 * 
 * @return 1 if it is integer, 0 otherwise
 */
int is_integer(double *xp)
{
    double x = *xp;
    double r = lround(x);
    if (fabs(r - x) < EPSILON)
    {
        *xp = r;
        return 1;
    } else
    {
        return 0;
    }
}

/**
 * @brief checks if elements in a node's solution vector are integers 
 * by iterating through the vector and calling is_integer.
 * 
 * @param p pointer to the node that we need to iterate through
 * 
 * @return 1 if all elements are integers, 0 otherwise
 */
int integer(node_t *p)
{
    for (int i = 0; i < p->n; i++)
    {
        if (!is_integer(&p->x[i]))
        {
            return 0;
        }
    }
    return 1;
}

/**
 * @brief updates current best solution and prunes the list of nodes.
 * 
 * the function compares a node's objective value with the current
 * best value and updates the solution if it is better. The function
 * also removes nodes from the list of unexplored nodes (h)
 * 
 * @param p pointer to the current node
 * @param h pointer to the head of the linked list of not yet explored nodes
 * @param zp pointer to the current best solution
 * @param x pointer to the solution vector that has the best value
 */
void bound(node_t* p, Node** h, double* zp, double* x){

    if(h == NULL || p == NULL){  return; }

    if(p->z > *zp)
    {
        *zp = p->z;
        
        memcpy(x, p->x, sizeof(double) * p->n);

        if (*h == NULL) { return; }

        Node *q, *prev, *next;
        q = *h;

        while (q->data->z < p->z)
        {
            q = q->next;

            if (q == NULL)      { return; }
            if (q->data == NULL)   { return; }
        }
        
        prev = q;
        q = q->next;

        while (q != NULL)
        {
            next = q->next;

            if (q->data->z < p->z)
            {
                prev->next = q->next;

            } else
            {
                prev = q;
            }
            q = next;

        }
    }
}

/**
 * @brief decides if the program should branch to a node.
 * 
 * @param q pointer to the node to evalaute
 * @param z current best solution.
 * 
 * @return 1 if branching possible, 0 otherwise.
 */
int branch(node_t* q, double z){
    double min, max;

    if (q->z < z) { return 0; }

    for (int h = 0; h < q->n; h++)
    {
        if (!is_integer(&q->x[h]))
        {
            if (q->min[h] == -INFINITY) { min = 0; }
            else                        { min = q->min[h]; }

            max = q->max[h];

            if (floor(q->x[h]) < min || ceil(q->x[h]) > max) { continue; }

            q->h = h;
            q->xh = q->x[h];

            return 1;
        }
    }
    return 0;
}

/**
 * @brief processes the succesor node in the branch and bound algorithm.
 * 
 * the succ function calls the extend function on the given node to extend
 * it with a new constraint. It then solves the subproblem, and either 
 * prunes, branches, or adds the node to the list of unexplored nodes.
 * After the node has been evalauted, all allocated memory is freed.
 * 
 * @param p Pointer to the parent node.
 * @param h Pointer to the head of the linked list of unexplored nodes.
 * @param m Number of constraints in the parent node.
 * @param n Number of decision variables in the parent node.
 * @param a Coefficient matrix of size m x n.
 * @param b RHS values of constraints (size m).
 * @param c Objective function coefficients (size n).
 * @param k Index of the decision variable for the new constraint.
 * @param ak Coefficient for the new constraint.
 * @param bk RHS value for the new constraint.
 * @param zp Pointer to the current best objective value (to be updated).
 * @param x Pointer to the solution vector corresponding to the best value.
 */
void succ (node_t* p, Node** h, int m, int n, double** a, double* b, double* c, int k, double ak, double bk, double* zp, double* x){
    node_t* q = extend(p, m, n, a, b, c, k, ak, bk);

    // if(q == NULL){ return; }

    q->z = simplex(q->m, q->n, q->a, q->b, q->c, q->x, 0);

    if(isfinite (q->z))
    {
        if(integer(q))
        {
            bound(q, h, zp, x);
        }else if (branch (q, *zp))
        {
            push(h, q);   
            return;     
        }
    }

    for (int i = 0; i < q->m + 1; i++)
    {
        free(q->a[i]);
        q->a[i] = NULL;
    }



    free(q->a);
    free(q->b);
    free(q->c);
    free(q->x);
    free(q->min);
    free(q->max);
    free(q);
    q = NULL;
}

/**
 * @brief solves the linear program problem through the branch and bound
 * algorithm
 * 
 * the branch and bound algorithm finds the optimal integer solution for
 * the linear program problem. The function intopt iteratively explores
 * and runes nodes while updating the best found solution.
 * 
 * @param m number of constraints
 * @param n number of decision variables
 * @param a coefficient matrix (m x n)
 * @param b b values of constraints
 * @param c objective function coefficents (n)
 * @param x pointer to the best solution vector
 * 
 * @return the optimal value of the objective function, or NAN if infeasible.
 */
double intopt(int m, int n, double** a, double* b, double* c, double* x)
{

    node_t* p   = initial_node(m,n,a,b,c);
    Node* h     = calloc(m, sizeof(node_t));

    h->data = p;
    double z = -INFINITY;

    p->z = simplex (p->m, p->n, p->a, p->b, p->c, p->x, 0);

    if(integer(p) || !isfinite(p->z))
    {
        z = p->z;
        if(integer(p)){
            memcpy(x, p->x, sizeof(double *) * p->n);
        }

        for (int i = 0; i < p->m + 1; i++)
        {
            free(p->a[i]);
        }

        free(p->a);
        free(p->b);
        free(p->c);
        free(p->x);
        free(p->min);
        free(p->max);
        free(p);
        p = NULL;

        return z;

    }

    branch(p,z);

    while(h != NULL){
        p = pop(&h);
        succ(p, &h, m, n, a, b, c, p->h, 1, floor(p->xh), &z, x);
        succ(p, &h, m, n, a, b, c, p->h, -1, -ceil(p->xh), &z, x);

        for (int i = 0; i < p->m + 1; i++)
        {
            free(p->a[i]);
        }

        free(p->a);
        free(p->b);
        free(p->c);
        free(p->x);
        free(p->min);
        free(p->max);
        free(p);
        p = NULL;
        
    }

    free(h);    

    if (z == -INFINITY)
    {
        return NAN;
    } else
    {
        return z;
    }
}