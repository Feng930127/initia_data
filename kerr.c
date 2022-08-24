static char help[] = "Einstein constraint equations\n\n";

#include <petsc.h>

#define EPSILON 1.0e-12

typedef struct {
  PetscReal u, v;       //u is conformal factor \psi,
} Field;                //v is stream function \Psi.

typedef struct {
  PetscReal  L,R,       //left and right side of the domain
             U,D,       //upper and down side of the domain
             a,M,
             coer, coethe;       //spin and mass of the black hole
} AppCtx;

/*------------------------------------------------------------------------
We define some functions here, rtuta is the quasi-isotropic radius for
the Kerr spacetime, similarily isotropic radius for the Scharzschild 
spacetime. 
*///----------------------------------------------------------------------
static PetscReal Delta_fun(PetscReal x, PetscReal y,
                           PetscReal a, PetscReal M) {
  PetscReal rtuta, Delta, r;
  
  r = PetscSqrtReal(x*x + y*y);
  rtuta = r + M + (M*M - a*a)/(4.0*r);
  Delta = a*a - 2.0*M*rtuta + rtuta*rtuta;

  if(PetscAbsReal(Delta) < EPSILON) return(0.);
  else return(Delta);
}


static PetscReal q_func(PetscReal x, PetscReal y,
                        PetscReal a, PetscReal M){
  PetscReal rtuta, Delta, r, the, Sigma,
            sinthe, costhe, sinthe2, costhe2, AA;

  r = PetscSqrtReal(x*x + y*y);
  the = PetscAsinReal(x/r); 
  rtuta = r + M + (M*M - a*a)/(4.0*r);
  Delta = a*a - 2.0*M*rtuta + rtuta*rtuta;
  sinthe  = PetscSinReal(the);
  costhe  = PetscCosReal(the);
  sinthe2 = PetscPowReal(sinthe,2.0);
  costhe2 = PetscPowReal(costhe,2.0);
  Delta = Delta_fun(x,y,a,M);
  Sigma = rtuta*rtuta + a*a*costhe2;
  AA    = PetscPowReal(a*a+rtuta*rtuta,2.0) - a*a*Delta*sinthe2;
  
  return (1./2.*PetscLogReal(Sigma*Sigma/AA));
}

static PetscReal CC_fun(PetscReal x, PetscReal y,
                        PetscReal a, PetscReal M) {
  PetscReal CC, r;
  
  r = PetscSqrtReal(x*x + y*y);
  CC = (1.0 - (M*M - a*a)/(r*r*r));

  if(PetscAbsReal(CC) < EPSILON) return(0.);
  else return(CC);
}

static PetscReal psi_exact(PetscReal x, PetscReal y,
                                PetscReal a, PetscReal M){
  PetscReal rtuta, sinthe, costhe, sinthe2,
            costhe2, Delta, Sigma, AA, BB, CC, r, the;

  r = PetscSqrtReal(x*x + y*y);
  the = PetscAsinReal(x/r);
  rtuta = r + M + (M*M - a*a)/(4.0*r);
  sinthe  = PetscSinReal(the);
  costhe  = PetscCosReal(the);
  sinthe2 = PetscPowReal(sinthe,2.0);
  costhe2 = PetscPowReal(costhe,2.0);
  Delta = Delta_fun(x,y,a,M);
  Sigma = rtuta*rtuta + a*a*costhe2;
  AA    = PetscPowReal(a*a+rtuta*rtuta,2.0) - a*a*Delta*sinthe2;
  CC    = CC_fun(x,y,a,M);
  BB    = 4.0*CC*rtuta*(a*a + rtuta*rtuta)
                - a*a*(-2.0*M*CC + 2.0*CC*rtuta)*sinthe2;

  return PetscPowReal(AA/(r*r*Sigma),1.0/4.0)
                + 0.0*(sinthe+costhe+costhe2+Delta+Sigma+CC+BB);
}

static PetscReal Omega_exact(PetscReal x, PetscReal y,
                                PetscReal a, PetscReal M){
  PetscReal rtuta, sinthe2, r, the, sinthe, AA, Omega;

  r = PetscSqrtReal(x*x + y*y);
  the = PetscAsinReal(x/r);
  rtuta = r + M + (M*M - a*a)/(4.0*r);
  sinthe  = PetscSinReal(the);
  sinthe2 = PetscPowReal(sinthe,2.0);
  AA    = PetscPowReal(a*a+rtuta*rtuta,2.0) - a*a*Delta_fun(x,y,a,M)*sinthe2;

  Omega = -2.0*a*M*rtuta/AA;
  return (Omega);
}

static PetscReal Psi_exact(PetscReal x, PetscReal y,
                                PetscReal a, PetscReal M){
  PetscReal rtuta, costhe, sinthe2,
            Sigma, r, the, costhe2;

  r = PetscSqrtReal(x*x + y*y);
  the = PetscAsinReal(x/r);
  rtuta = r + M + (M*M - a*a)/(4.0*r);
  costhe  = PetscCosReal(the);
  sinthe2 = PetscPowReal(PetscSinReal(the),2.0);
  costhe2 = PetscPowReal(PetscCosReal(the),2.0);
  Sigma = rtuta*rtuta + a*a*costhe2;

  return 2.0*M*a*(costhe2*costhe-3.0*costhe)-2.0*M*a*a*a*costhe*sinthe2*sinthe2/Sigma;
}

static PetscReal dNdr(PetscReal x, PetscReal y,
                                PetscReal a, PetscReal M){
  PetscReal rtuta, sinthe2, dNdr, r, the,
            costhe2, Delta, Sigma, AA, BB, CC;
 
  r = PetscSqrtReal(x*x + y*y);
  the = PetscAsinReal(x/r);
  rtuta = r + M + (M*M - a*a)/(4.0*r);
  sinthe2 = PetscPowReal(PetscSinReal(the),2.0);
  costhe2 = PetscPowReal(PetscCosReal(the),2.0);
  Delta = Delta_fun(x,y,a,M);
  Sigma = rtuta*rtuta + a*a*costhe2;
  AA    = PetscPowReal(a*a+rtuta*rtuta,2.0)-a*a*Delta*sinthe2;
  CC    = CC_fun(x,y,a,M);
  BB    = 4.0*CC*rtuta*(a*a + rtuta*rtuta)
                - a*a*(-2.0*M*CC + 2.0*CC*rtuta)*sinthe2;

  dNdr = - (AA*AA*CC*(Delta*rtuta + (M - rtuta)*Sigma)
             + a*a*rtuta*(-AA*AA*CC+4.0*AA*CC*M*M*(rtuta*rtuta-Sigma)
             + 2.0*BB*M*M*rtuta*Sigma)*sinthe2)/(AA*Sigma*(AA*Delta
             - a*a*(AA - 4.0*M*M*rtuta*rtuta)*sinthe2));

  if(PetscAbsReal(dNdr) < EPSILON) return(0.);
  else return(dNdr);
}

static PetscReal dNdtheta(PetscReal x, PetscReal y,
                                PetscReal a, PetscReal M){
  PetscReal rtuta, sinthe, costhe, sinthe2, r, the,
            costhe2, Delta, Sigma, AA, dNdtheta;

  r = PetscSqrtReal(x*x + y*y);
  the = PetscAsinReal(x/r);
  rtuta = r + M + (M*M - a*a)/(4.0*r);
  sinthe  = PetscSinReal(the);
  costhe  = PetscCosReal(the);
  sinthe2 = PetscPowReal(PetscSinReal(the),2.0);
  costhe2 = PetscPowReal(PetscCosReal(the),2.0);
  Delta = Delta_fun(x,y,a,M);
  Sigma = rtuta*rtuta + a*a*costhe2;
  AA    = PetscPowReal(a*a+rtuta*rtuta,2.0) - a*a*Delta*sinthe2;

  dNdtheta = -a*a*costhe*sinthe/Sigma + a*a*Delta*costhe*sinthe/AA;

  if(PetscAbsReal(dNdtheta) < EPSILON) return(0.);
  else return(dNdtheta);
}
static PetscReal dNdro(PetscReal x, PetscReal y,
                        PetscReal a, PetscReal M){
  PetscReal the, ans, r;
  
  r = PetscSqrtReal(x*x + y*y);
  the = PetscAsinReal(x/r);
  ans = dNdr(x,y,a,M)*PetscSinReal(the) 
        + dNdtheta(x,y,a,M)*PetscCosReal(the);
  return(ans);
}
static PetscReal dNdz(PetscReal x, PetscReal y,
                        PetscReal a, PetscReal M){
  PetscReal the, ans, r;
  
  r = PetscSqrtReal(x*x + y*y);
  the = PetscAsinReal(x/r);
  ans = dNdr(x,y,a,M)*PetscCosReal(the) 
        - dNdtheta(x,y,a,M)*PetscSinReal(the)/r;
  return(ans);
}
static PetscReal NN(PetscReal x, PetscReal y,
                        PetscReal a, PetscReal M){
  PetscReal r, the, sinthe2, costhe2, Sigma,
            AA, rtuta, Delta, ans;  
  r = PetscSqrtReal(x*x + y*y);
  the = PetscAsinReal(x/r);
  rtuta = r + M + (M*M - a*a)/(4.0*r);
  sinthe2 = PetscPowReal(PetscSinReal(the),2.0);
  costhe2 = PetscPowReal(PetscCosReal(the),2.0);
  Delta = Delta_fun(x,y,a,M);
  Sigma = rtuta*rtuta + a*a*costhe2;
  AA    = PetscPowReal(a*a+rtuta*rtuta,2.0)-a*a*Delta*sinthe2;

  ans= PetscSqrtReal((4.0*a*a*M*M*rtuta*rtuta*sinthe2*sinthe2)/(Sigma*Sigma) 
		+ (AA*sinthe2*(Delta - a*a*sinthe2))/(Sigma*Sigma))
		/PetscSqrtReal((AA*sinthe2)/Sigma);
 
  return (ans);
}

static PetscReal dqdrho(PetscReal x, PetscReal y,
                        PetscReal a, PetscReal M){
  PetscReal rtuta, sinthe, costhe, sinthe2, r, the,
            Delta, Sigma, AA, BB, CC, costhe2, ans, dqdr,dqdtheta;

  r = PetscSqrtReal(x*x + y*y);
  the = PetscAsinReal(x/r);
  rtuta = r + M + (M*M - a*a)/(4.0*r);
  sinthe  = PetscSinReal(the);
  costhe  = PetscCosReal(the);
  sinthe2 = PetscPowReal(PetscSinReal(the),2.0);
  costhe2 = PetscPowReal(PetscCosReal(the),2.0);
  Delta = Delta_fun(x,y,a,M);
  Sigma = rtuta*rtuta + a*a*costhe2;
  AA    = PetscPowReal(a*a+rtuta*rtuta,2.0)-a*a*Delta*sinthe2;
  CC    = (1.0 - (M*M - a*a)/(4*r*r));
  BB    = 4.0*CC*rtuta*(a*a + rtuta*rtuta)
                - a*a*(-2.0*M*CC + 2.0*CC*rtuta)*sinthe2;
  dqdr = AA*(4.0*CC*rtuta*Sigma/AA-BB*Sigma*Sigma/AA/AA)/(2.0*Sigma*Sigma);
  dqdtheta = AA*(-4.0*a*a*Sigma*costhe*sinthe/AA+2.0*a*a*Delta*Sigma*Sigma*costhe*sinthe/(AA*AA))/(2.0*Sigma*Sigma);
  ans = dqdr*sinthe + dqdtheta*costhe/r;

  return(ans); 

}
static PetscReal dqdz(PetscReal x, PetscReal y,
                        PetscReal a, PetscReal M){
  PetscReal rtuta, sinthe, costhe, sinthe2, r, the,
            Delta, Sigma, AA, BB, CC, costhe2, ans, dqdr,dqdtheta;

  r = PetscSqrtReal(x*x + y*y);
  the = PetscAsinReal(x/r);
  rtuta = r + M + (M*M - a*a)/(4.0*r);
  sinthe  = PetscSinReal(the);
  costhe  = PetscCosReal(the);
  sinthe2 = PetscPowReal(PetscSinReal(the),2.0);
  costhe2 = PetscPowReal(PetscCosReal(the),2.0);
  Delta = Delta_fun(x,y,a,M);
  Sigma = rtuta*rtuta + a*a*costhe2;
  AA    = PetscPowReal(a*a+rtuta*rtuta,2.0)-a*a*Delta*sinthe2;
  CC    = (1.0 - (M*M - a*a)/(4*r*r));
  BB    = 4.0*CC*rtuta*(a*a + rtuta*rtuta)
                - a*a*(-2.0*M*CC + 2.0*CC*rtuta)*sinthe2;
  dqdr = AA*(4.0*CC*rtuta*Sigma/AA-BB*Sigma*Sigma/AA/AA)/(2.0*Sigma*Sigma);
  dqdtheta = AA*(-4.0*a*a*Sigma*costhe*sinthe/AA+2.0*a*a*Delta*Sigma*Sigma*costhe*sinthe/(AA*AA))/(2.0*Sigma*Sigma);
  ans = dqdr*costhe - dqdtheta*sinthe/r;

  return(ans); 

}
static PetscReal qro2z2(PetscReal x, PetscReal y,
                                PetscReal a, PetscReal M){
  PetscReal rtuta, sinthe, costhe, sinthe2, r, the,
            Delta, Sigma, AA, BB, CC, costhe2;

  r = PetscSqrtReal(x*x + y*y);
  the = PetscAsinReal(x/r);
  rtuta = r + M + (M*M - a*a)/(4.0*r);
  sinthe  = PetscSinReal(the);
  costhe  = PetscCosReal(the);
  sinthe2 = PetscPowReal(PetscSinReal(the),2.0);
  costhe2 = PetscPowReal(PetscCosReal(the),2.0);
  Delta = Delta_fun(x,y,a,M);
  Sigma = rtuta*rtuta + a*a*costhe2;
  AA    = PetscPowReal(a*a+rtuta*rtuta,2.0)-a*a*Delta*sinthe2;
  CC    = (1.0 - (M*M - a*a)/(4*r*r));
  BB    = 4.0*CC*rtuta*(a*a + rtuta*rtuta)
                - a*a*(-2.0*M*CC + 2.0*CC*rtuta)*sinthe2;

 return BB*(-BB*Sigma*Sigma/(AA*AA) + 4.0*CC*rtuta*Sigma/AA)/(2.0*Sigma*Sigma)
        - (2.0*AA*CC*rtuta*((4.0*CC*rtuta*Sigma)/AA-BB*Sigma*Sigma/(AA*AA)))/(Sigma*Sigma*Sigma)
        + (AA*(-BB*Sigma*Sigma/(AA*AA) + 4.0*CC*rtuta*Sigma/AA))/(2*Sigma*Sigma*x)
        + (1.0/(2.0*Sigma*Sigma))*AA*(-((8.0*BB*CC*Sigma*rtuta)/(AA*AA))
                + (8.0*CC*CC*rtuta*rtuta)/AA + (4.0*CC*CC*Sigma)/AA
                + (2.0*BB*BB*Sigma*Sigma)/(AA*AA*AA)
                + (2.0*(-a*a + M*M)*rtuta*Sigma)/(AA*r*r*r)
                - (Sigma*Sigma*(8.0*CC*CC*rtuta*rtuta + 4.0*CC*CC*(a*a + rtuta*rtuta)
                + (2.0*(-a*a + M*M)*rtuta*(a*a + rtuta*rtuta))/(r*r*r)
                - a*a*(2.0*CC*CC - (M*(-a*a + M*M))/(r*r*r)
                + ((-a*a + M*M)*rtuta)/(r*r*r))*sinthe2))/(AA*AA))
        + (1.0/(r*r))*((2.0*a*a*AA*costhe*sinthe*(-((4.0*a*a*Sigma*costhe*sinthe)/AA)
                + (2.0*a*a*Delta*Sigma*Sigma*costhe*sinthe)/(AA*AA)))/(Sigma*Sigma*Sigma)
                - (a*a*Delta*costhe*sinthe*(-((4.0*a*a*Sigma*costhe*sinthe)/AA)
                + (2.0*a*a*Delta*Sigma*Sigma*costhe*sinthe)/(AA*AA)))/(Sigma*Sigma)
                + (1.0/(2.0*Sigma*Sigma))*AA*(-((4.0*a*a*Sigma*costhe2)/AA)
                + (2.0*a*a*Delta*Sigma*Sigma*costhe2)/(AA*AA)
                + (4.0*a*a*Sigma*sinthe2)/AA - (2.0*a*a*Delta*Sigma*Sigma*sinthe2)/(AA*AA)
                + (8.0*a*a*a*a*costhe2*sinthe2)/AA - (16.0*a*a*a*a*Delta*Sigma*costhe2*sinthe2)/(AA*AA)
                + (8.0*a*a*a*a*Delta*Delta*Sigma*Sigma*costhe2*sinthe2)/(AA*AA*AA)));
}
//-------------------------------------------------------------------------



//Barrier functions for robin boundary conditions-----------------
static PetscReal barrier(PetscReal x, PetscReal y) {
  PetscReal barrier, a, r2;
  
  a = PetscPowReal(10.0,10.0);
  r2 = PetscSqrtReal((x - 0.001)*(x - 0.001) + (y)*(y)) - 1.5;
  barrier = - 1.5/PETSC_PI*PetscAtanReal(r2*a) + 1.5/2;
 
  if(PetscAbsReal(barrier) < 0.1) return(0.);
  else return(barrier);
}

static PetscReal barrier_rup(PetscReal x, PetscReal y, 
			   PetscReal r0, PetscReal epsilon) {
  PetscReal barrier_rup, a, r1, r2;

  a = PetscPowReal(10.0, 10.0);
  r1 = PetscSqrtReal((x - 0.001)*(x - 0.001) + (y)*(y)) 
		- (r0 - epsilon);
  r2 = PetscSqrtReal((x - 0.001)*(x - 0.001) + (y)*(y)) 
		- (r0 + epsilon);
  barrier_rup = 1.0/PETSC_PI*(PetscAtanReal(a*r1)
		- PetscAtanReal(a*r2)); 
  
  return(barrier_rup);
}

static PetscReal barrier_thetaup(PetscReal x, PetscReal y,
			       PetscReal theta0, PetscReal epsilon) {
  PetscReal barrier_thetaup, a, theta1, theta2;
  
  a = PetscPowReal(10.0, 10.0);
  theta1 = PetscAtanReal((y)/x) - (theta0 - epsilon);
  theta2 = PetscAtanReal((y)/x) - (theta0 + epsilon);
  barrier_thetaup = 1.0/PETSC_PI*(PetscAtanReal(a*theta1)
                - PetscAtanReal(a*theta2));
  
  return (barrier_thetaup);
}
//-----------------------------------------------------------------------


//Define some functions-------------------------------------------------- 
extern PetscErrorCode InitialState(DM, Vec, AppCtx*);
extern PetscErrorCode ExactSolution(DM, Vec, AppCtx*);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo*, Field**, Field**, AppCtx*);
extern PetscErrorCode FormJacobianLocal(DMDALocalInfo*, Field**, Mat, Mat, AppCtx*);
//------------------------------------------------------------------------


//STARTMAIN
int main(int argc,char **argv) {
  PetscErrorCode ierr;
  DM            da;
  SNES          snes;
  AppCtx        user;
  Vec           x, Yexact, err, err1, y;
  Mat           M1;
  PetscReal     hx, hy, errnorm, errnorm1, errnorm2;
  DMDALocalInfo info;

  ierr = PetscInitialize(&argc,&argv,NULL,help); if (ierr) return ierr;
  user.a  = 0.1;
  user.M  = 1.0;
  user.L  = 0.001;
  user.R  = 10.0;
  user.U  = 10.0;
  user.D  = -10.0; 
  user.coer = 0.52;
  user.coethe = 0.52; 
  
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD, "Kerr_",
                                "options for Kerr", ""); CHKERRQ(ierr);
  ierr = PetscOptionsReal("-a","The spin of the black hole",
           "jacobian_neumann_constraints_kerr_Psi.c",user.a,&user.a,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-M","The mass of the black hole",
           "jacobian_neumann_constraints_kerr_Psi.c",user.M,&user.M,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-L","The left side of the domain",
           "jacobian_neumann_constraints_kerr_Psi.c",user.L,&user.L,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-R","The right side of the domain",
           "jacobian_neumann_constraints_kerr_Psi.c",user.R,&user.R,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-U","The upper side of the domain",
           "jacobian_neumann_constraints_kerr_Psi.c",user.U,&user.U,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-D","The down side of the domain",
           "jacobian_neumann_constraints_kerr_Psi.c",user.D,&user.D,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  ierr = DMDACreate2d(PETSC_COMM_WORLD,
               DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
               DMDA_STENCIL_STAR,  // for 5-point stencil
               41, 41, PETSC_DECIDE,PETSC_DECIDE,
               2, 1,              // degrees of freedom, stencil width
               NULL, NULL, &da); CHKERRQ(ierr);
  ierr = DMSetFromOptions(da); CHKERRQ(ierr);
  ierr = DMSetUp(da); CHKERRQ(ierr);
  ierr = DMDASetFieldName(da, 0, "u"); CHKERRQ(ierr);
  ierr = DMDASetFieldName(da, 1, "v"); CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(da, &info); CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(da, user.L, user.R,
                                user.D, user.U, -1.0, -1.0); CHKERRQ(ierr);
  ierr = DMSetApplicationContext(da,&user); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
           "Solving Einstein constraint equations on %d x %d grid\n", 
			info.mx,info.my); CHKERRQ(ierr);

  ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(da,&x); CHKERRQ(ierr);
  ierr = VecDuplicate(x,&Yexact); CHKERRQ(ierr);
  ierr = VecDuplicate(x,&err); CHKERRQ(ierr);
  ierr = VecDuplicate(x,&err1); CHKERRQ(ierr);
  ierr = VecDuplicate(x,&y); CHKERRQ(ierr);
  ierr = InitialState(da,x,&user); CHKERRQ(ierr);
  ierr = ExactSolution(da,Yexact,&user); CHKERRQ(ierr);
  ierr = DMCreateMatrix(da,&M1); CHKERRQ(ierr);

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetDM(snes,da); CHKERRQ(ierr);
  ierr = DMDASNESSetFunctionLocal(da,INSERT_VALUES,
             (DMDASNESFunction)FormFunctionLocal,&user); CHKERRQ(ierr);
  //ierr = DMDASNESSetJacobianLocal(da,
    //         (DMDASNESJacobian)FormJacobianLocal,&user); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
  
  ierr = SNESSolve(snes,NULL,x); CHKERRQ(ierr);
  ierr = SNESComputeFunction(snes, x, y); CHKERRQ(ierr);
 
  hx = (user.R - user.L) / (PetscReal)(info.mx - 1);
  hy = (user.U - user.D) / (PetscReal)(info.my - 1);
  ierr = VecWAXPY(err,-1.0,x,Yexact); CHKERRQ(ierr);
  ierr = VecWAXPY(err1,-1.0,x,Yexact); CHKERRQ(ierr);
  ierr = VecScale(err1, hx*hy); CHKERRQ(ierr);
  ierr = VecNorm(err, NORM_INFINITY, &errnorm); CHKERRQ(ierr);
  ierr = VecNorm(err1, NORM_1, &errnorm1); CHKERRQ(ierr);

  ierr = VecNorm(y, NORM_INFINITY, &errnorm2); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,
           "L-Inf(err)=%g, L1(err)=%g, L-Inf(res)=%g\n", errnorm,errnorm1,errnorm2); CHKERRQ(ierr);

  VecDestroy(&x);    VecDestroy(&Yexact); 
  VecDestroy(&err);   VecDestroy(&err1);
  SNESDestroy(&snes);    DMDestroy(&da);
  return PetscFinalize();
}
//ENDMAIN

//Define the initial state----------------------------------------------
PetscErrorCode InitialState(DM da, Vec Y, AppCtx* user) {
  PetscErrorCode   ierr;
  DMDALocalInfo    info;
  PetscInt         i,j;
  PetscReal        hx, hy, x, y, a = user->a, M = user->M;
  DMDACoor2d       **aC;
  Field            **aY;

  ierr = VecSet(Y,0.0); CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);
  ierr = DMDAGetCoordinateArray(da,&aC); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,Y,&aY); CHKERRQ(ierr);
  hx = (user->R - user->L) / (PetscReal)(info.mx - 1);
  hy = (user->U - user->D) / (PetscReal)(info.my - 1);

  for (j = info.ys; j < info.ys+info.ym; j++) {
      y = j * hy + user->D;
      for (i = info.xs; i < info.xs+info.xm; i++) {
          x = i * hx + user->L;
          if (i==info.mx-1||j==0||j==info.my-1) {
                aY[j][i].u = psi_exact(x, y, a, M);
                aY[j][i].v = Omega_exact(x, y, a, M);
             }
             else{
                aY[j][i].u = 1.0 + 0.001*(1.5 - barrier(x, y))
                                 + 0.5*barrier(x, y);
		        aY[j][i].v = 0.1*barrier(x, y);

             }
      }
  }
  ierr = DMDAVecRestoreArray(da,Y,&aY); CHKERRQ(ierr);
  ierr = DMDARestoreCoordinateArray(da,&aC); CHKERRQ(ierr);
  return 0;
}

//Exact Solutions--------------------------------------------------------
PetscErrorCode ExactSolution(DM da, Vec Y, AppCtx* user) {
  PetscErrorCode   ierr;
  DMDALocalInfo    info;
  PetscInt         i,j;
  PetscReal        x, y, hx, hy,
                   a = user->a, M = user->M;
  DMDACoor2d       **aC;
  Field            **aY;

  ierr = VecSet(Y,0.0); CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);
  ierr = DMDAGetCoordinateArray(da,&aC); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,Y,&aY); CHKERRQ(ierr);
  hx = (user->R - user->L) / (PetscReal)(info.mx - 1);
  hy = (user->U - user->D) / (PetscReal)(info.my - 1);

  for (j = info.ys; j < info.ys+info.ym; j++) {
      y = hy * j + user->D;
      for (i = info.xs; i < info.xs+info.xm; i++) {
          x = hx * i + user->L;
          aY[j][i].u = psi_exact(x, y, a, M);
          aY[j][i].v = Omega_exact(x, y, a, M);
      }
  }
  ierr = DMDAVecRestoreArray(da,Y,&aY); CHKERRQ(ierr);
  ierr = DMDARestoreCoordinateArray(da,&aC); CHKERRQ(ierr);
  return 0;
}

//index function for y > 0 part------------------------------------------
PetscErrorCode indexfunc_up(DMDALocalInfo *info, PetscReal theta,
                        PetscReal r1, PetscReal r2, PetscReal *value, 
  			Field **aY, AppCtx *user) {
  PetscInt   i, j, mx = info->mx, my = info->my, k = 0, l = 0;
  PetscReal  hx = (user->R - user->L) / (PetscReal)(mx - 1),
             hy = (user->U - user->D) / (PetscReal)(my - 1),
             x, y, coer = user->coer, coethe = user->coethe, 
             barrier1, barrier2;

  for (j = info->ys; j < info->ys + info->ym; j++) {
      y = hy * j + user->D;
      for (i = info->xs; i < info->xs + info->xm; i++) {
          x = hx * i + user->L; 
          barrier1 = barrier_rup(x, y, r1, coer*hx)
                        *barrier_thetaup(x, y, theta, coethe*hx);
 	      barrier2 = barrier_rup(x, y, r2, coer*hx)
                        *barrier_thetaup(x, y, theta, coethe*hx);
	  if (barrier1 > 0.9) {
                value[0] = aY[j][i].u;
		k = k + 1;
          }
	  else if (barrier2 > 0.9) {
		value[1] = aY[j][i].u;
		l = l + 1;
	  }
	  else if (k == 1 && l == 1) {
		goto here;
	  }
      }
  }
  here:
  return 0.;
}


//Define the residual function--------------------------------------------
PetscErrorCode FormFunctionLocal(DMDALocalInfo *info, Field **aY, 
				Field **aG, AppCtx *user) {
  PetscErrorCode ierr;
  PetscInt   i, j, mx = info->mx, my = info->my, n, k;
  PetscReal  value[2],uxx, uyy, u1, u2, u3, vxx, vyy, v1, v2, v3, x, y;
  PetscReal  hx = (user->R - user->L) / (PetscReal)(mx -1),
             hy = (user->U - user->D) / (PetscReal)(my -1),
             a = user->a, M = user->M, vxl, vxr, vyl, vyr,
             uxl, uxr, uyl, uyr, theta, barrier_rup1,
             barrier_thetaup1, part, coer = user->coer;

  for (j = info->ys; j < info->ys + info->ym; j++) {
      y = hy * j + user->D;
      for (i = info->xs; i < info->xs + info->xm; i++) {
	  x = hx * i + user->L; theta = 0.0;
	  barrier_rup1 = barrier_rup(x, y, 1.5 - coer*hx, coer*hx);
	  if (barrier_rup1 > 0.9) { 
	      n = (PetscInt)(PETSC_PI/(coer*hx));
	      for (k = 0; k < n + 1; k++) {
	   	  theta = (2.0 * k)*coer*hx - PETSC_PI/2.0;
	  	  barrier_thetaup1 = barrier_thetaup(x, y, theta, coer*hx);
		  if (barrier_thetaup1 > 0.6) {
		      ierr = indexfunc_up(info, theta, 1.5 + 3.0*coer*hx,
				 1.5 + coer*hx, value, aY, user); CHKERRQ(ierr);
		      part = 1.0/(value[1]*PetscExpReal(q_func(x,y,a,M))*1.5);
		      aY[j][i].u = value[0] + coer*hx*part*(x*dqdrho(x,y,a,M)
				+ (y)*dqdz(x,y,a,M) + 3.0);
		      break;
 	  	   }
	       }
	   }
      }
  }
	
  for (j = info->ys; j < info->ys + info->ym; j++) {
      y = hy * j + user->U;
      for (i = info->xs; i < info->xs + info->xm; i++) {
          x = hx * i + user->L;
          if (i==0||i==mx-1||j==0||j==my-1) {
              aG[j][i].u =0.0;// (aY[j][i].u - psi_exact(x, y, a, M));
              aG[j][i].v = 0.0;//aY[j][i].v - Omega_exact(x, y, a, M));
          }
          else {
	          uxl = (i==0) ? aY[j][i+1].u : aY[j][i-1].u;
              uxr = (i==mx-1) ? aY[j][i-1].u : aY[j][i+1].u;
              uyl = (j==0) ? aY[j+1][i].u : aY[j-1][i].u;
              uyr = (j==my-1) ? aY[j-1][i].u : aY[j+1][i].u;
              vxl = (i==0) ? aY[j][i+1].v : aY[j][i-1].v;
              vxr = (i==mx-1) ? aY[j][i-1].v : aY[j][i+1].v;
              vyl = (j==0) ? aY[j+1][i].v : aY[j-1][i].v;
              vyr = (j==my-1) ? aY[j-1][i].v : aY[j+1][i].v;

              uxx = (uxl - 2.0 * aY[j][i].u + uxr)/(hx*hx);
              uyy = (uyl - 2.0 * aY[j][i].u + uyr)/(hy*hy);
              u1 = (aY[j][i+1].u - uxl) / (2.0*hx*x);
              u2 = 1.0/4.0*qro2z2(x,y,a,M)*aY[j][i].u;
              u3 = 1.0/16.0*PetscPowReal(NN(x,y,a,M), -2.0)*x*x
                     * (PetscPowReal(1.0/(2.0*hx)*(aY[j][i+1].v - vxl),2.0)
                         + PetscPowReal(1.0/(2.0*hy)*(aY[j+1][i].v - aY[j-1][i].v),2.0))
                     * PetscPowReal(aY[j][i].u,5.0);

              vxx = (vxl - 2.0 * aY[j][i].v + vxr)/(hx*hx);
              vyy = (vyl - 2.0 * aY[j][i].v + vyr)/(hy*hy);
               v1 = 3.0/x*1.0/(2.0*hx)*(aY[j][i+1].v - vxl);
              v2 = -dNdro(x,y,a,M)*1.0/(2.0*hx)*(aY[j][i+1].v - vxl)
                        - dNdz(x,y,a,M)*1.0/(2.0*hy)*(aY[j+1][i].v - aY[j-1][i].v);
              v3 = 6.0/(aY[j][i].u)*1.0/(2.0*hx)*(aY[j][i+1].u - uxl)
                                *1.0/(2.0*hx)*(aY[j][i+1].v - vxl)
                   + 6.0/(aY[j][i].u)*1.0/(2.0*hy)*(aY[j+1][i].u - aY[j-1][i].u)
                                *1.0/(2.0*hy)*(aY[j+1][i].v - aY[j-1][i].v);
          //    v1 = -3.0/x*1.0/(2.0*hx)*(aY[j][i+1].v - vxl);
          //    v2 = -dNdro(x,y)*1.0/(2.0*hx)*(aY[j][i+1].v - vxl)
            //            - dNdz(x,y)*1.0/(2.0*hy)*(aY[j+1][i].v - aY[j-1][i].v);
              //v3 = -6.0/(aY[j][i].u)*1.0/(2.0*hx)*(aY[j][i+1].u - uxl)
                //                *1.0/(2.0*hx)*(aY[j][i+1].v - vxl)
                  // - 6.0/(aY[j][i].u)*1.0/(2.0*hy)*(aY[j+1][i].u - aY[j-1][i].u)
                    //            *1.0/(2.0*hy)*(aY[j+1][i].v - aY[j-1][i].v);
	          aG[j][i].u = (uxx + uyy + 1.0*(u1 + u2 + u3))*(1.5 - barrier(x, y))
                                + 0.0*barrier(x, y);
              aG[j][i].v = (vxx + vyy + 1.0*(v1 + v2 + v3))*(1.5 - barrier(x,y))
                                + 0.0*barrier(x, y);
	 
          }
      }
    }
    return 0;
}

PetscErrorCode FormJacobianLocal(DMDALocalInfo *info, Field **aY,
                                    Mat J, Mat P, AppCtx *user) {
  PetscErrorCode ierr;
  PetscInt    i, j, mx = info->mx, my = info->my;
  PetscReal   v[5],L=user->L, R=user->R, U=user->U, D=user->D;
  MatStencil  col[5],row;
  PetscReal   hr = (user->R - user->L) / (PetscReal)(mx -1),
              htheta = (user->D - user->U) / (PetscReal)(my -1),
              a = user->a, M = user->M, vil, vir, vjl, vjr,
              uil, uir, ujl, ujr, x, y;
  
  for (j = info->ys; j < info->ys+info->ym; j++) {
      row.j = j;  y = htheta * j + user->U;
      for (i = info->xs; i < info->xs+info->xm; i++) {
          row.i = i;  x = hr * i + user->L;
	  if (i==0||i==mx-1||j==0||j==my-1) {
	      row.c = 0;  col[0].c = 0; 
              col[0].i = i; col[0].j = j; v[0] = 1.0;
	      ierr = MatSetValuesStencil(P,1,&row,1,col,v,INSERT_VALUES); CHKERRQ(ierr);
	      row.c = 1;  col[0].c = 1; 
              col[0].i = i; col[0].j = j; v[0] = 1.0;
              ierr = MatSetValuesStencil(P,1,&row,1,col,v,INSERT_VALUES); CHKERRQ(ierr);
          }
	  else {
  	  uil = (i==1) ? psi_exact(L, y, a, M) : aY[j][i-1].u;
          uir = (i==mx-2) ? psi_exact(R, y, a, M) : aY[j][i+1].u;
          ujl = (j==1) ? psi_exact(x, U, a, M) : aY[j-1][i].u;
          ujr = (j==my-2) ? psi_exact(x, D, a, M) : aY[j+1][i].u;
          vil = (i==1) ? Psi_exact(L, y, a, M) : aY[j][i-1].v;
          vir = (i==mx-2) ? Psi_exact(R, y, a, M) : aY[j][i+1].v;
          vjl = (j==1) ? Psi_exact(x, U, a, M) : aY[j-1][i].v;
          vjr = (j==my-2) ? Psi_exact(x, D, a, M) : aY[j+1][i].v;
	  //psi equation
	  row.c = 0;  col[0].c = 0;  col[1].c = 0;
	  col[2].c = 0;	col[3].c = 0; col[4].c = 0;
          col[0].i = i; col[0].j = j; 
	  v[0] = - 2.0*(htheta/hr + hr/htheta/(x*x)) + 1.0/4.0*qro2z2(x,y,a,M)*hr*htheta
	  	    + 1.0/16.0*PetscPowReal(x*PetscSinReal(y),-4.0)
                       * (PetscPowReal(1.0/(2.0)*(vir-vil),2.0)*htheta/hr
                           + 1.0/(x*x)*PetscPowReal(1.0/(2.0)*(vjr-vjl),2.0)*hr/htheta)
                       * PetscPowReal(aY[j][i].u,-8.0)*(-7.0);
	  col[1].i = i-1; col[1].j = j; 
	  v[1] = (i>1) ? 1.0*htheta/hr - 1.0/(x)*htheta : 0.0; 
	  col[2].i = i+1; col[2].j = j; 
	  v[2] = (i<mx-2) ? 1.0*htheta/hr + 1.0/(x)*htheta : 0.0;
   	  col[3].i = i; col[3].j = j-1; 
	  v[3] = (j>1) ? 1.0/(x*x)*hr/htheta - 1.0/(2.0*x*x)*hr/PetscTanReal(y) : 0.0;
   	  col[4].i = i; col[4].j = j+1; 
	  v[4] = (j<my-2) ? 1.0/(x*x)*hr/htheta + 1.0/(2.0*x*x)*hr/PetscTanReal(y) : 0.0;
   	  ierr = MatSetValuesStencil(P,1,&row,5,col,v,INSERT_VALUES); CHKERRQ(ierr);
	  
	  row.c = 0;  col[0].c = 1;  col[1].c = 1;
          col[2].c = 1; col[3].c = 1; col[4].c = 1;
	  col[0].i = i; col[0].j = j;
 	  v[0]= 0.0;
	  col[1].i = i-1; col[1].j = j;
	  v[1] = (i>1) ? 1.0/16.0*PetscPowReal(x*PetscSinReal(y),-4.0)
                   * (-1.0*PetscPowReal(1.0/2.0*(vir-vil),1.0)*htheta/hr)
                   	* PetscPowReal(aY[j][i].u,-7.0) : 0.0;
	  col[2].i = i+1; col[2].j = j;
          v[2] = (i<mx-2) ? 1.0/16.0*PetscPowReal(x*PetscSinReal(y),-4.0)
                   * (1.0*PetscPowReal(1.0/2.0*(vir-vil),1.0)*htheta/hr)
                        * PetscPowReal(aY[j][i].u,-7.0) : 0.0;
 	  col[3].i = i; col[3].j = j-1;
   	  v[3] = (j>1) ? 1.0/16.0*PetscPowReal(x*PetscSinReal(y),-4.0)
                   * (- 1.0/(x*x)*PetscPowReal(1.0/(2.0)*(vjr-vjl),1.0)*hr/htheta)
                   *PetscPowReal(aY[j][i].u,-7.0) : 0.0;
	  col[4].i = i; col[4].j = j+1;
          v[4] = (j<my-2) ? 1.0/16.0*PetscPowReal(x*PetscSinReal(y),-4.0)
                   * (1.0/(x*x)*PetscPowReal(1.0/(2.0)*(vjr-vjl),1.0)*hr/htheta)
                   *PetscPowReal(aY[j][i].u,-7.0) : 0.0;
	  ierr = MatSetValuesStencil(P,1,&row,5,col,v,INSERT_VALUES); CHKERRQ(ierr);


 	  //Y equation
	  row.c = 1;  col[0].c = 1;  col[1].c = 1;
          col[2].c = 1; col[3].c = 1; col[4].c = 1;
          col[0].i = i; col[0].j = j; 
	  v[0] = - 2.0*(htheta/hr + hr/htheta/(x*x));
 	  col[1].i = i-1; col[1].j = j; 
	  v[1] = (i>1) ? 1.0*htheta/hr + 2.0/(2.0*x)*htheta - dNdr(x,y,a,M)/2.0*htheta 
			+ 6.0/(2.0*aY[j][i].u)*htheta/(2.0*hr)*(uir-uil) : 0.0;
	  col[2].i = i+1; col[2].j = j;
	  v[2] = (i<mx-2) ? 1.0*htheta/hr - 2.0/(2.0*x)*htheta + dNdr(x,y,a,M)/2.0*htheta
			- 6.0/(2.0*aY[j][i].u)*htheta/(2.0*hr)*(uir-uil) : 0.0;
	  col[3].i = i; col[3].j = j-1;
  	  v[3] = (j>1) ? 1.0/(x*x)*hr/htheta + 3.0/(2.0*x*x*PetscTanReal(y))*hr
			-dNdtheta(x,y,a,M)/(2.0*x*x)*hr
			+ 6.0/(2.0*aY[j][i].u*x*x)*hr/(2.0*htheta)*(ujr-ujl) : 0.0;
	  col[4].i = i; col[4].j = j+1;
	  v[4] = (j<my-2) ? 1.0/(x*x)*hr/htheta - 3.0/(2.0*x*x*PetscTanReal(y))*hr
                	+ dNdtheta(x,y,a,M)/(2.0*x*x)*hr 
			- 6.0/(2.0*aY[j][i].u*x*x)*hr/(2.0*htheta)*(ujr-ujl) : 0.0;
	  ierr = MatSetValuesStencil(P,1,&row,5,col,v,INSERT_VALUES); CHKERRQ(ierr);

	  row.c = 1;  col[0].c = 0;  col[1].c = 0;
          col[2].c = 0; col[3].c = 0; col[4].c = 0;
          col[0].i = i; col[0].j = j;
	  v[0] = 6.0/(aY[j][i].u*aY[j][i].u)*1.0/(2.0)*(uir-uil)
                        *1.0/(2.0*hr)*(vir-vil)*htheta
                + 6.0/(aY[j][i].u*aY[j][i].u)*1.0/(2.0*x*x)*(ujr-ujl)
                        *1.0/(2.0*htheta)*(vjr-vjl)*hr;
	  col[1].i = i-1; col[1].j = j;
	  v[1] = (i>1) ? - 6.0/(aY[j][i].u)*1.0/(2.0)*(-1.0)
                                *1.0/(2.0*hr)*(vir - vil)*htheta : 0.0;
	  col[2].i = i+1; col[2].j = j;
	  v[2] = (i<mx-2) ? - 6.0/(aY[j][i].u)*1.0/(2.0)*(1.0)
                                *1.0/(2.0*hr)*(vir - vil)*htheta : 0.0;
	  col[3].i = i; col[3].j = j-1;
	  v[3] = (j>1) ? - 6.0/(aY[j][i].u)*1.0/(2.0*x*x)*(-1.0)
                                *1.0/(2.0*htheta)*(vjr - vjl)*hr : 0.0;
	  col[3].i = i; col[4].j = j+1;
 	  v[4] = (j<my-2) ? - 6.0/(aY[j][i].u)*1.0/(2.0*x*x)*(1.0)
                                *1.0/(2.0*htheta)*(vjr - vjl)*hr : 0.0;
	  ierr = MatSetValuesStencil(P,1,&row,5,col,v,INSERT_VALUES); CHKERRQ(ierr);
	  }
        }
    }
  ierr = MatAssemblyBegin(P,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(P,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  //ierr = MatView(P,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  if (J != P) {
      ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  }
  return 0;
}
