
//MPI Types
void Create_rot_datatype(MPI_Datatype *Rottype);
void Create_side_datatype(MPI_Datatype *Sidetype);
void Create_back_datatype(MPI_Datatype *Backtype);
void Create_pot_datatype(MPI_Datatype *Pottype);
void Create_parms_datatype(MPI_Datatype *Parmstype);

//Send-Receive
struct s_mc_parms *send_parms(int iproc, int nprocs, MPI_Datatype Parmstype, struct s_mc_parms *parms, MPI_Status astatus);//,double *tmin);
void send_struct(int *nback, int iproc, int nprocs, int *nat, int *ntypes, MPI_Status astatus);
struct s_polymer *send_pol(int iproc, int nprocs, int nback, MPI_Datatype Backtype, MPI_Datatype Sidetype,  MPI_Datatype Rottype, struct s_polymer *startp, MPI_Status astatus, int npol, int shell, int nosidechains);
struct s_potential *send_pot(int nat, int ntypes, int noangpot, int nodihpot, int hb, int iproc, int nprocs, MPI_Datatype Pottype, struct s_potential *u, MPI_Status astatus);
void send_double_matrix(int length1, int length2, int iproc, double **m, int source);
void send_int_matrix(int length1, int length2, int iproc, int **m, int source);

//Exchange
void ExchangePol(struct s_polymer *polymer, struct s_polymer *replica, struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *u, int iproc, int ntemp, int even, int *ex_count, int *ex_acc, MPI_Datatype Backtype, MPI_Datatype Sidetype, MPI_Datatype Rottype, MPI_Status astatus);

//Potential-related
#ifdef OPTIMIZEPOT
struct s_optimizepot_input *Allo_op_input(int nrestr);
struct s_optimizepot_input *send_op_input(int nrestr,struct s_optimizepot_input *in);
#endif
