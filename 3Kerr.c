static char help[] = "Einstein constraint equations\n\n";

#include <petsc.h>

#define EPSILON 1.0e-12

typedef struct {
  PetscReal u, v;       //u is conformal factor \psi,
} Field;                //v is stream function \Psi.

typedef struct {
  PetscReal  L,R,       //left and right side of the domain
             U,D,       //upper and down side of the domain
             a,M,       //spin and mass of the black hole
             coer, coethe,
             rhori;       
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
  PetscReal rtuta, Delta, r, the, Sigma,ans,
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
  ans =  PetscPowReal(AA/(r*r*Sigma),1.0/4.0);
  
  if(PetscAbsReal(ans) < EPSILON) return(0.);
  else return(ans);
}
static PetscReal CC_fun(PetscReal x, PetscReal y,
                         PetscReal a, PetscReal M) {
   PetscReal CC, r;
   
   r = PetscSqrtReal(x*x + y*y);
   CC = (1.0 - (M*M - a*a)/(4.0*r*r));

  if(PetscAbsReal(CC) < EPSILON) return(0.);
  else return(CC);
 }
static PetscReal psi_exact(PetscReal x, PetscReal y,
                                  PetscReal a, PetscReal M){
  PetscReal rtuta, sinthe, costhe, sinthe2, ans,
            costhe2, Delta, Sigma, AA, r, the;

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
  ans =  PetscPowReal(AA/(r*r*Sigma),1.0/4.0);
  
  if(PetscAbsReal(ans) < EPSILON) return(0.);
  else return(ans);
}

static PetscReal Omega_exact(PetscReal x, PetscReal y,
                                PetscReal a, PetscReal M){
  PetscReal rtuta, sinthe2, r, the, sinthe, AA, Omega;

  r = PetscSqrtReal(x*x + y*y);
  the = PetscAsinReal(x/r);
  rtuta = r + M + (M*M - a*a)/(4.0*r);
  sinthe  = PetscSinReal(the);
  sinthe2 = PetscPowReal(sinthe,2.0);
  AA    = PetscPowReal(a*a+rtuta*rtuta,2.0) 
                - a*a*Delta_fun(x,y,a,M)*sinthe2;
  Omega = - 2.0*a*M*rtuta/AA;
  
  if(PetscAbsReal(Omega) < EPSILON) return(0.);
  else return(Omega);
}

static PetscReal Psi_exact(PetscReal x, PetscReal y,
                                PetscReal a, PetscReal M){
  PetscReal rtuta, costhe, sinthe2,
            Sigma, r, the, costhe2, ans;

  r = PetscSqrtReal(x*x + y*y);
  the = PetscAsinReal(x/r);
  rtuta = r + M + (M*M - a*a)/(4.0*r);
  costhe  = PetscCosReal(the);
  sinthe2 = PetscPowReal(PetscSinReal(the),2.0);
  costhe2 = PetscPowReal(PetscCosReal(the),2.0);
  Sigma = rtuta*rtuta + a*a*costhe2;
  ans =  2.0*M*a*(costhe2*costhe - 3.0*costhe) 
                - 2.0*M*a*a*a*costhe*sinthe2*sinthe2/Sigma;

  if(PetscAbsReal(ans) < EPSILON) return(0.);
  else return(ans);
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
  dNdtheta = 1.0/PetscTanReal(the) - 1.0/PetscSqrtReal(AA*sinthe2/Sigma)
 		* ((2.0*AA*costhe*sinthe)/Sigma + (2.0*a*a*AA*costhe*sinthe2*sinthe)/(Sigma*Sigma)
                - (2.0*a*a*Delta*costhe*sinthe2*sinthe)/Sigma)/(2.0*PetscSqrtReal((AA*sinthe2)/Sigma));

  if(PetscAbsReal(dNdtheta) < EPSILON) return(0.);
  else return(dNdtheta);
}
static PetscReal dNdro(PetscReal x, PetscReal y,
                        PetscReal a, PetscReal M){
  PetscReal the, ans, r;
  
  r = PetscSqrtReal(x*x + y*y);
  the = PetscAsinReal(x/r);
  ans = dNdr(x,y,a,M)*PetscSinReal(the) 
        + dNdtheta(x,y,a,M)*PetscCosReal(the)/r;
//printf("x=%g, y=%g, r=%g,the=%g,ans=%g\n", x, y, r, the, ans);
  if(PetscAbsReal(ans) < EPSILON) return(0.);
  else return(ans);
 
}
static PetscReal dNdz(PetscReal x, PetscReal y,
                        PetscReal a, PetscReal M){
  PetscReal the, ans, r;
  
  r = PetscSqrtReal(x*x + y*y);
  the = PetscAsinReal(x/r);
  ans = dNdr(x,y,a,M)*PetscCosReal(the) 
        - dNdtheta(x,y,a,M)*PetscSinReal(the)/r;
 // printf("x=%g, y=%g, r=%g,the=%g,ans=%g\n", x, y, r, the, ans);
  if(PetscAbsReal(ans) < EPSILON) return(0.);
  else return(ans);
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
  if(PetscAbsReal(ans) < EPSILON) return(0.);
  else return(ans);

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
  dqdtheta = AA*(- 4.0*a*a*Sigma*costhe*sinthe/AA
            + 2.0*a*a*Delta*Sigma*Sigma*costhe*sinthe/(AA*AA))/(2.0*Sigma*Sigma);
  ans = dqdr*sinthe + dqdtheta*costhe/r;
   
  if(PetscAbsReal(ans) < EPSILON) return(0.);
  else return(ans);

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
  dqdtheta = AA*(- 4.0*a*a*Sigma*costhe*sinthe/AA
            + 2.0*a*a*Delta*Sigma*Sigma*costhe*sinthe/(AA*AA))/(2.0*Sigma*Sigma);
  ans = dqdr*costhe - dqdtheta*sinthe/r;
   
  if(PetscAbsReal(ans) < EPSILON) return(0.);
  else return(ans);
}
static PetscReal qro2z2(PetscReal x, PetscReal y,
                                PetscReal a, PetscReal M){
  PetscReal ans, rtuta, sinthe, costhe, sinthe2, r, the,
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

 ans =  BB*(-BB*Sigma*Sigma/(AA*AA) + 4.0*CC*rtuta*Sigma/AA)/(2.0*Sigma*Sigma)
        - (2.0*AA*CC*rtuta*((4.0*CC*rtuta*Sigma)/AA-BB*Sigma*Sigma/(AA*AA)))/(Sigma*Sigma*Sigma)
        + (AA*(-BB*Sigma*Sigma/(AA*AA) + 4.0*CC*rtuta*Sigma/AA))/(2*Sigma*Sigma*r)
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
//printf("x=%g, y=%g, r=%g,the=%g,rtuta=%g,ans=%g\n", x, y, r, the, rtuta,ans);  
  if(PetscAbsReal(ans) < EPSILON) return(0.);
  else return(ans);
}
//-------------------------------------------------------------------------



//Barrier functions for robin boundary conditions-----------------
static PetscReal barrier(PetscReal x, PetscReal y, 
                            PetscReal rhori) {
  PetscReal barrier, A, r1;
  
  A = PetscPowReal(10.0,22.0);
  r1 = PetscSqrtReal((x - 0.0001)*(x - 0.0001) + (y)*(y)) - (rhori);
  barrier = - 1.0/PETSC_PI*PetscAtanReal(r1*A) + 1.0/2;
  if(PetscAbsReal(barrier) < 0.1) return(0.);
  else return(barrier);
}

static PetscReal barrier_rup(PetscReal x, PetscReal y, 
			   PetscReal r0, PetscReal epsilon) {
  PetscReal barrier_rup, A, r1, r2;

  A = PetscPowReal(10.0, 22.0);
  r1 = PetscSqrtReal((x - 0.0001)*(x - 0.0001) + (y)*(y)) 
		- (r0 - epsilon);
  r2 = PetscSqrtReal((x - 0.0001)*(x - 0.0001) + (y)*(y)) 
		- (r0 + epsilon);
  barrier_rup = 1.0/PETSC_PI*(PetscAtanReal(A*r1)
		- PetscAtanReal(A*r2)); 
  
  return(barrier_rup);
}

static PetscReal barrier_thetaup(PetscReal x, PetscReal y,
			       PetscReal theta0, PetscReal epsilon) {
  PetscReal barrier_thetaup, a, theta1, theta2;
  
  a = PetscPowReal(10.0, 22.0);
  theta1 = PetscAtanReal((y)/(x-0.0001)) - (theta0 - epsilon);
  theta2 = PetscAtanReal((y)/(x-0.0001)) - (theta0 + epsilon);
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
extern PetscErrorCode err_func(DMDALocalInfo*, DM, Vec, Vec, AppCtx*);
//------------------------------------------------------------------------


//STARTMAIN
int main(int argc,char **argv) {
  PetscErrorCode ierr;
  DM            da;
  SNES          snes;
  AppCtx        user;
  Vec           x, Yexact, err, y;
  Mat           M1;
  PetscReal     hx, hy, errnorm, errnorm1, errnorm2;
  DMDALocalInfo info;

  ierr = PetscInitialize(&argc,&argv,NULL,help); if (ierr) return ierr;
  user.a  = 0.3;
  user.M  = 3.0;
  user.L  = 0.0001;
  user.R  = 20.0;
  user.U  = 10.0;
  user.D  = -10.0; 
  user.coer = 0.50;
  user.coethe = 0.50; 
  user.rhori = PetscSqrtReal(user.M*user.M - user.a*user.a)/2.0;
 
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
               5, 5, PETSC_DECIDE,PETSC_DECIDE,
               2, 1,              // degrees of freedom, stencil width
               NULL, NULL, &da); CHKERRQ(ierr);
  ierr = DMSetFromOptions(da); CHKERRQ(ierr);
  ierr = DMSetUp(da); CHKERRQ(ierr);
  ierr = DMDASetFieldName(da, 0, "psi"); CHKERRQ(ierr);
  ierr = DMDASetFieldName(da, 1, "Omega"); CHKERRQ(ierr);
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
  ierr = err_func(&info, da, x, Yexact, &user); CHKERRQ(ierr); 

  hx = (user.R - user.L) / (PetscReal)(info.mx - 1);
  hy = (user.U - user.D) / (PetscReal)(info.my - 1);
  ierr = VecWAXPY(err,-1.0,x,Yexact); CHKERRQ(ierr);
  ierr = VecNorm(err, NORM_INFINITY, &errnorm); CHKERRQ(ierr);
  ierr = VecNorm(err, NORM_1, &errnorm1); CHKERRQ(ierr);

  ierr = VecNorm(y, NORM_INFINITY, &errnorm2); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,
           "L-Inf(err)=%g, L1(err)=%g, L-Inf(res)=%g\n", errnorm,errnorm1*hx*hy,errnorm2/hx/hy); CHKERRQ(ierr);

  VecDestroy(&x);    VecDestroy(&Yexact); 
  VecDestroy(&err);  
  SNESDestroy(&snes);    DMDestroy(&da);
  return PetscFinalize();
}
//ENDMAIN

//Define the initial state----------------------------------------------
PetscErrorCode InitialState(DM da, Vec Y, AppCtx* user) {
  PetscErrorCode   ierr;
  DMDALocalInfo    info;
  PetscInt         i,j;
  PetscReal        hx, hy, x, y, a = user->a, M = user->M,
                   rhori = user->rhori, coer = user->coer;
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
          if (i == info.mx-1 || j == 0 || j == info.my-1) {
                aY[j][i].u = psi_exact(x, y, a, M);
                aY[j][i].v = Omega_exact(x, y, a, M);
             }
             else{
                aY[j][i].u = 1.0 + 1.00434*barrier(x, y, rhori);
		        aY[j][i].v = Omega_exact(x,y,a,M)*(1.0 - barrier(x,y,rhori))
                             - 0.00835427*barrier(x,y,rhori);

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
                   a = user->a, M = user->M,
                   rhori = user->rhori;
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
          aY[j][i].u = psi_exact(x, y, a, M)*(1.0 - barrier(x, y, rhori)) 
                            + 2.00434*barrier(x, y, rhori);
          aY[j][i].v = Omega_exact(x, y, a, M)*(1.0 - barrier(x, y, rhori)) 
                            - 0.00835427*barrier(x, y, rhori);
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
	  if (barrier1 > 0.8) {
                value[0] = aY[j][i].u;
		k = k + 1;
          }
	  else if (barrier2 > 0.8) {
		value[1] = aY[j][i].u;
		l = l + 1;
	  }
	  else if (k == 1 && l == 1) {
		goto here;
	  }
      }
  }
  here:
  if (l == 0 || k == 0){
    PetscPrintf(PETSC_COMM_WORLD,"l =%g, k =%g\n", l, k); 
  }

  return 0.;
}


//Define the residual function--------------------------------------------
PetscErrorCode FormFunctionLocal(DMDALocalInfo *info, Field **aY, 
				Field **aG, AppCtx *user) {
  PetscErrorCode ierr;
  PetscInt   i, j, mx = info->mx, my = info->my;
  PetscReal  value[2], uxx, uyy, u1, u2, u3, vxx, vyy, v1, v2, v3, x, y;
  PetscReal  hx = (user->R - user->L) / (PetscReal)(mx - 1),
             hy = (user->U - user->D) / (PetscReal)(my - 1),
             a = user->a, M = user->M, vxl, vxr, vyl, vyr,
             uxl, uxr, uyl, uyr, theta, barrier_rup1,
             barrier_thetaup1, part, coer = user->coer,
             rhori = user->rhori;
  //for (j = info->ys; j < info->ys + info->ym; j++) {
  //    y = hy * j + user->D;
  //    for (i = info->xs; i < info->xs + info->xm; i++) {
	//  x = hx * i + user->L; theta = 0.0;
//	  barrier_rup1 = barrier_rup(x, y, rhori - 1.0*coer*hx, coer*hx);
 // 	  if (barrier_rup1 > 0.5) { 
  //        theta = PetscAtanReal(y/(x - 0.0001));
//		  ierr = indexfunc_up(info, theta, rhori + 3.0*coer*hx,
//				 rhori + 1.0*coer*hx, value, aY, user); CHKERRQ(ierr);
//		  part = value[1]/(2.0*(rhori + 1.0*coer*hx));
//		  aY[j][i].u = (4.0*value[1] - value[0] + 4.0*coer*hx*part)/3.0;
//          break;
//	   }
//     }
//  }

  for (j = info->ys; j < info->ys + info->ym; j++) {
      y = hy * j + user->D;
      for (i = info->xs; i < info->xs + info->xm; i++) {
	  x = hx * i + user->L;
	  barrier_rup1 = barrier_rup(x, y, rhori - 0.0*coer*hx, coer*hx);
      //if (i==0){
      //     aG[j][i].u = -3.0*aY[j][i].u + 4.0*aY[j][i+1].u -aY[j][i+2].u;
      //     aG[j][i].v = -3.0*aY[j][i].v + 4.0*aY[j][i+1].v -aY[j][i+2].v;
      //}
      
	  if (barrier_rup1 > 0.6) { 
	      aG[j][i].v = aY[j][i].v - Omega_exact(x, y, a, M);
	   	  theta = PetscAtanReal(y/(x - 0.0001));
	      ierr = indexfunc_up(info, theta, rhori + 4.0*coer*hx,
	   	            rhori + 2.0*coer*hx, value, aY, user); CHKERRQ(ierr);
	      part = -3.0*aY[j][i].u + 4.0*value[1] - value[0];
	      aG[j][i].u = part/(4.0*coer*hx) + aY[j][i].u/2.0/(rhori);
          //aG[j][i].u = aY[j][i].u - psi_exact(x, y, a, M);
       }
       else if (i==mx-1||j==0||j==my-1&barrier_rup1<=0.6){
           aG[j][i].u = aY[j][i].u - psi_exact(x, y, a, M);
           aG[j][i].v = aY[j][i].v - Omega_exact(x, y, a, M);
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

           uxx = (uxl - 2.0 * aY[j][i].u + uxr)/(hx)*hy;
           uyy = (uyl - 2.0 * aY[j][i].u + uyr)/(hy)*hx;
           u1 = (aY[j][i+1].u - uxl) / (2.0*x)*hy;
           u2 = 1.0/4.0*qro2z2(x,y,a,M)*aY[j][i].u*hx*hy;
           u3 = 1.0/16.0*PetscPowReal(NN(x,y,a,M), -2.0)*x*x
                * (PetscPowReal(1.0/(2.0)*(aY[j][i+1].v - vxl),2.0)*hy/hx
                    + PetscPowReal(1.0/(2.0)*(aY[j+1][i].v - aY[j-1][i].v),2.0)*hx/hy)
                * PetscPowReal(aY[j][i].u, 5.0);

           vxx = (vxl - 2.0 * aY[j][i].v + vxr)/(hx)*hy;
           vyy = (vyl - 2.0 * aY[j][i].v + vyr)/(hy)*hx;
           v1 = 3.0/x*1.0/(2.0)*(aY[j][i+1].v - vxl)*hy;
           v2 = - dNdro(x,y,a,M)*1.0/(2.0)*(aY[j][i+1].v - vxl)*hy
                - dNdz(x,y,a,M)*1.0/(2.0)*(aY[j+1][i].v - aY[j-1][i].v)*hx;
           v3 =  6.0/(aY[j][i].u)*1.0/(2.0)*(aY[j][i+1].u - uxl)
                        *1.0/(2.0)*(aY[j][i+1].v - vxl)*hy/hx
                + 6.0/(aY[j][i].u)*1.0/(2.0)*(aY[j+1][i].u - aY[j-1][i].u)
                        *1.0/(2.0)*(aY[j+1][i].v - aY[j-1][i].v)*hx/hy;
	          
           aG[j][i].u = (uxx + uyy + u1 + u2 + u3)*(1.0 - barrier(x, y, rhori))
                            + 0.0*barrier(x, y, rhori);
           aG[j][i].v = 0.0*(vxx + vyy + v1 + v2 + v3)*(1.0 - barrier(x, y, rhori))
                            + 0.0*barrier(x, y, rhori);
       } 
    }
  }
	
    return 0;
}


PetscErrorCode err_func(DMDALocalInfo *info, DM da, Vec Y, 
                    Vec Yexact, AppCtx *user) {
  
  PetscErrorCode   ierr;
  PetscInt         i, j;
  PetscReal        x, y, hx, hy, a = user->a, M = user->M, 
                   err2 = 0.0, err1 = 0.0, mx = info->mx, 
                   my = info->my, rhori = user->rhori;
  DMDACoor2d       **aC;
  Field            **aY, **aYexact;

  ierr = DMDAGetCoordinateArray(da,&aC); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,Y,&aY); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, Yexact, &aYexact); CHKERRQ(ierr);
  hx = (user->R - user->L) / (PetscReal)(mx - 1);
  hy = (user->U - user->D) / (PetscReal)(my - 1);
 
  for (j = info->ys; j < info->ys + info->ym; j++) {
      y = hy * j + user->D;
      for (i = info->xs; i < info->xs + info->xm; i++) {
	      x = hx * i + user->L;  
          err1 = err1 + PetscAbsReal(aY[j][i].u - aYexact[j][i].u)
                            * (1.0 - barrier(x, y, rhori))*hx*hy;
          err2 = err2 + PetscAbsReal(aY[j][i].v - aYexact[j][i].v)
                            * (1.0 - barrier(x, y, rhori))*hx*hy;
      }
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "err1=%g\n", err1); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,Y,&aY); CHKERRQ(ierr);
  ierr = DMDARestoreCoordinateArray(da,&aC); CHKERRQ(ierr);
  return 0;
 
}

