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
             sig1,sig2,
             coer,coethe;       //spin and mass of the black hole
} AppCtx;

static PetscReal Power(PetscReal x, PetscReal y) {
  PetscReal value;

  value = PetscPowReal(x, y);
  if (PetscAbsReal(value) < EPSILON) return (0.);
  else return (value);
}

static PetscReal Sqrt(PetscReal x) {
  PetscReal value;

  value = PetscSqrtReal(x);
  if (PetscAbsReal(value) < EPSILON) return (0.);
  else return (value);
}

static PetscReal dOmegadrho(PetscReal x, PetscReal y,
                        PetscReal M, PetscReal a) {
  PetscReal value;
  value = (a*M*x*(Power(a,2) - Power(M,2) + 4*(Power(x,2) + Power(y,2)))*
     (3*Power(a,8) - 4*Power(a,6)*(3*Power(M,2) + 8*Power(x,2) - 4*Power(y,2) + 
          12*M*Sqrt(Power(x,2) + Power(y,2))) + 
       3*(Power(M,8) + 16*Power(M,7)*Sqrt(Power(x,2) + Power(y,2)) + 
          112*Power(M,6)*(Power(x,2) + Power(y,2)) + 
          448*Power(M,5)*Power(Power(x,2) + Power(y,2),1.5) + 
          1120*Power(M,4)*Power(Power(x,2) + Power(y,2),2) + 
          1792*Power(M,3)*Power(Power(x,2) + Power(y,2),2.5) + 
          1792*Power(M,2)*Power(Power(x,2) + Power(y,2),3) + 
          1024*M*Power(Power(x,2) + Power(y,2),3.5) + 
          256*Power(Power(x,2) + Power(y,2),4)) + 
       2*Power(a,4)*(9*Power(M,4) + 72*Power(M,3)*Sqrt(Power(x,2) + Power(y,2)) + 
          32*M*Sqrt(Power(x,2) + Power(y,2))*(7*Power(x,2) + 3*Power(y,2)) + 
          8*Power(M,2)*(25*Power(x,2) + 19*Power(y,2)) - 
          16*(-5*Power(x,4) + 2*Power(x,2)*Power(y,2) + 7*Power(y,4))) - 
       4*Power(a,2)*(3*Power(M,6) + 36*Power(M,5)*Sqrt(Power(x,2) + Power(y,2)) + 
          448*M*Power(Power(x,2) + Power(y,2),2.5) + 
          64*Power(Power(x,2) + Power(y,2),2)*(2*Power(x,2) + 3*Power(y,2)) + 
          64*Power(M,3)*Sqrt(Power(x,2) + Power(y,2))*
           (7*Power(x,2) + 6*Power(y,2)) + 
          4*Power(M,4)*(44*Power(x,2) + 41*Power(y,2)) + 
          48*Power(M,2)*(13*Power(x,4) + 24*Power(x,2)*Power(y,2) + 11*Power(y,4)))))
    /(512.*Power(Power(x,2) + Power(y,2),3.5)*
     Power(-0.0625*(Power(a,2)*Power(x,2)*
           Power(Power(a,2) - Power(M,2) + 4*(Power(x,2) + Power(y,2)),2))/
         Power(Power(x,2) + Power(y,2),2) + 
       Power(Power(a,2) + Power(M + 
           (-Power(a,2) + Power(M,2))/(4.*Sqrt(Power(x,2) + Power(y,2))) + 
           Sqrt(Power(x,2) + Power(y,2)),2),2),2));
  if (PetscAbsReal(value) < EPSILON) return (0.);
  else return (value);
}

static PetscReal dOmegadz(PetscReal x, PetscReal y,
                        PetscReal M, PetscReal a) {
  PetscReal value;
  value = (a*M*y*(Power(a,2) - Power(M,2) + 4*(Power(x,2) + Power(y,2)))*
     (16*Power(a,2)*Power(x,2)*Power(Power(a,2) - Power(M,2) + 
          4*(Power(x,2) + Power(y,2)),2) - 
       Power(16*Power(a,2)*(Power(x,2) + Power(y,2)) + 
         Power(-Power(a,2) + Power(M,2) + 4*M*Sqrt(Power(x,2) + Power(y,2)) + 
           4*(Power(x,2) + Power(y,2)),2),2) - 
       4*(Power(a,2) - Power(M,2) - 4*M*Sqrt(Power(x,2) + Power(y,2)) - 
          4*(Power(x,2) + Power(y,2)))*
        (8*Power(a,2)*Power(x,2)*(Power(a,2) - Power(M,2) - 
             4*(Power(x,2) + Power(y,2))) + 
          8*Power(a,2)*Power(x,2)*(Power(a,2) - Power(M,2) + 
             4*(Power(x,2) + Power(y,2))) - 
          (Power(a,2) - Power(M,2) - 4*M*Sqrt(Power(x,2) + Power(y,2)) - 
             4*(Power(x,2) + Power(y,2)))*
           (16*Power(a,2)*(Power(x,2) + Power(y,2)) + 
             Power(-Power(a,2) + Power(M,2) + 4*M*Sqrt(Power(x,2) + Power(y,2)) + 
               4*(Power(x,2) + Power(y,2)),2)))))/
   (512.*Power(Power(x,2) + Power(y,2),3.5)*
     Power(-0.0625*(Power(a,2)*Power(x,2)*
           Power(Power(a,2) - Power(M,2) + 4*(Power(x,2) + Power(y,2)),2))/
         Power(Power(x,2) + Power(y,2),2) + 
       Power(Power(a,2) + Power(M + 
           (-Power(a,2) + Power(M,2))/(4.*Sqrt(Power(x,2) + Power(y,2))) + 
           Sqrt(Power(x,2) + Power(y,2)),2),2),2));
 if (PetscAbsReal(value) < EPSILON) return (0.);
 else return (value);
}
static PetscReal LapOmega12(PetscReal x, PetscReal y, 
                        PetscReal M, PetscReal a) {
  PetscReal value;
  value = (256*a*M*Power(x,3)*Sqrt(Power(x,2) + Power(y,2))*
     (9*Power(a,18) - 3*Power(a,16)*
        (27*Power(M,2) + 32*(Power(x,2) - Power(y,2)) + 
          88*M*Sqrt(Power(x,2) + Power(y,2))) + 
       3*M*(-3*Power(M,17) - 88*Power(M,16)*Sqrt(Power(x,2) + Power(y,2)) - 
          1184*Power(M,15)*(Power(x,2) + Power(y,2)) - 
          9600*Power(M,14)*Power(Power(x,2) + Power(y,2),1.5) - 
          51520*Power(M,13)*Power(Power(x,2) + Power(y,2),2) - 
          186368*Power(M,12)*Power(Power(x,2) + Power(y,2),2.5) - 
          419328*Power(M,11)*Power(Power(x,2) + Power(y,2),3) - 
          292864*Power(M,10)*Power(Power(x,2) + Power(y,2),3.5) + 
          1830400*Power(M,9)*Power(Power(x,2) + Power(y,2),4) + 
          8785920*Power(M,8)*Power(Power(x,2) + Power(y,2),4.5) + 
          22257664*Power(M,7)*Power(Power(x,2) + Power(y,2),5) + 
          38764544*Power(M,6)*Power(Power(x,2) + Power(y,2),5.5) + 
          49201152*Power(M,5)*Power(Power(x,2) + Power(y,2),6) + 
          45875200*Power(M,4)*Power(Power(x,2) + Power(y,2),6.5) + 
          30801920*Power(M,3)*Power(Power(x,2) + Power(y,2),7) + 
          14155776*Power(M,2)*Power(Power(x,2) + Power(y,2),7.5) + 
          3997696*M*Power(Power(x,2) + Power(y,2),8) + 
          524288*Power(Power(x,2) + Power(y,2),8.5)) - 
       4*Power(a,14)*(-81*Power(M,4) - 
          528*Power(M,3)*Sqrt(Power(x,2) + Power(y,2)) + 
          32*M*Sqrt(Power(x,2) + Power(y,2))*(-19*Power(x,2) + 13*Power(y,2)) - 
          48*Power(M,2)*(22*Power(x,2) + 15*Power(y,2)) + 
          48*(-Power(x,4) + 6*Power(x,2)*Power(y,2) + 3*Power(y,4))) + 
       Power(a,2)*(81*Power(M,16) + 2112*Power(M,15)*Sqrt(Power(x,2) + Power(y,2)) + 
          589824*Power(Power(x,2) + Power(y,2),8) + 
          2621440*M*Power(Power(x,2) + Power(y,2),6.5)*(2*Power(x,2) + Power(y,2)) + 
          786432*Power(M,2)*Power(Power(x,2) + Power(y,2),6)*
           (15*Power(x,2) + 2*Power(y,2)) + 
          18432*Power(M,9)*Power(Power(x,2) + Power(y,2),2.5)*
           (20*Power(x,2) + 11*Power(y,2)) - 
          1179648*Power(M,3)*Power(Power(x,2) + Power(y,2),5.5)*
           (-Power(x,2) + 12*Power(y,2)) + 
          192*Power(M,14)*(130*Power(x,2) + 129*Power(y,2)) - 
          98304*Power(M,5)*Power(Power(x,2) + Power(y,2),4.5)*
           (797*Power(x,2) + 747*Power(y,2)) - 
          49152*Power(M,7)*Power(Power(x,2) + Power(y,2),3.5)*
           (1070*Power(x,2) + 913*Power(y,2)) + 
          128*Power(M,13)*Sqrt(Power(x,2) + Power(y,2))*
           (1369*Power(x,2) + 1337*Power(y,2)) - 
          49152*Power(M,6)*Power(Power(x,2) + Power(y,2),4)*
           (1688*Power(x,2) + 1463*Power(y,2)) - 
          16384*Power(M,4)*Power(Power(x,2) + Power(y,2),5)*
           (2315*Power(x,2) + 2827*Power(y,2)) - 
          4608*Power(M,8)*Power(Power(x,2) + Power(y,2),3)*
           (3889*Power(x,2) + 3377*Power(y,2)) + 
          1024*Power(M,10)*Power(Power(x,2) + Power(y,2),2)*
           (3901*Power(x,2) + 3476*Power(y,2)) + 
          576*Power(M,12)*(1385*Power(x,4) + 2706*Power(x,2)*Power(y,2) + 
             1321*Power(y,4)) + 1536*Power(M,11)*Sqrt(Power(x,2) + Power(y,2))*
           (1539*Power(x,4) + 2963*Power(x,2)*Power(y,2) + 1424*Power(y,4))) + 
       4*Power(a,12)*(-189*Power(M,6) - 
          1848*Power(M,5)*Sqrt(Power(x,2) + Power(y,2)) - 
          336*Power(M,4)*(20*Power(x,2) + 17*Power(y,2)) - 
          96*Power(M,3)*Sqrt(Power(x,2) + Power(y,2))*
           (113*Power(x,2) + 49*Power(y,2)) + 
          128*M*Sqrt(Power(x,2) + Power(y,2))*
           (-Power(x,4) + 31*Power(x,2)*Power(y,2) + 20*Power(y,4)) + 
          48*Power(M,2)*(-135*Power(x,4) - 38*Power(x,2)*Power(y,2) + 
             77*Power(y,4)) - 384*(-Power(x,6) - Power(x,4)*Power(y,2) + 
             9*Power(x,2)*Power(y,4) + 9*Power(y,6))) + 
       2*Power(a,10)*(567*Power(M,8) + 
          7392*Power(M,7)*Sqrt(Power(x,2) + Power(y,2)) + 
          672*Power(M,6)*(58*Power(x,2) + 53*Power(y,2)) + 
          960*Power(M,5)*Sqrt(Power(x,2) + Power(y,2))*
           (109*Power(x,2) + 77*Power(y,2)) - 
          768*Power(Power(x,2) + Power(y,2),2)*
           (7*Power(x,4) - 18*Power(x,2)*Power(y,2) + 55*Power(y,4)) + 
          480*Power(M,4)*(293*Power(x,4) + 378*Power(x,2)*Power(y,2) + 
             93*Power(y,4)) - 256*Power(M,3)*Sqrt(Power(x,2) + Power(y,2))*
           (-253*Power(x,4) - 29*Power(x,2)*Power(y,2) + 176*Power(y,4)) - 
          512*Power(M,2)*(59*Power(x,6) + 206*Power(x,4)*Power(y,2) + 
             187*Power(x,2)*Power(y,4) + 40*Power(y,6)) + 
          1024*M*Sqrt(Power(x,2) + Power(y,2))*
           (-28*Power(x,6) + 3*Power(x,4)*Power(y,2) + 78*Power(x,2)*Power(y,4) + 
             47*Power(y,6))) - 4*Power(a,6)*
        (-189*Power(M,12) - 3696*Power(M,11)*Sqrt(Power(x,2) + Power(y,2)) - 
          336*Power(M,10)*(94*Power(x,2) + 91*Power(y,2)) - 
          480*Power(M,9)*Sqrt(Power(x,2) + Power(y,2))*
           (319*Power(x,2) + 287*Power(y,2)) + 
          12288*Power(Power(x,2) + Power(y,2),4)*
           (-Power(x,4) + 6*Power(x,2)*Power(y,2) + 3*Power(y,4)) + 
          8192*M*Power(Power(x,2) + Power(y,2),3.5)*
           (Power(x,4) + 20*Power(x,2)*Power(y,2) + 7*Power(y,4)) - 
          4096*Power(M,2)*Power(Power(x,2) + Power(y,2),3)*
           (-116*Power(x,4) - 127*Power(x,2)*Power(y,2) + 13*Power(y,4)) + 
          12288*Power(M,3)*Power(Power(x,2) + Power(y,2),2.5)*
           (124*Power(x,4) + 137*Power(x,2)*Power(y,2) + 17*Power(y,4)) + 
          768*Power(M,4)*Power(Power(x,2) + Power(y,2),2)*
           (2443*Power(x,4) + 3302*Power(x,2)*Power(y,2) + 979*Power(y,4)) - 
          240*Power(M,8)*(1871*Power(x,4) + 3350*Power(x,2)*Power(y,2) + 
             1483*Power(y,4)) - 256*Power(M,7)*Sqrt(Power(x,2) + Power(y,2))*
           (2933*Power(x,4) + 4765*Power(x,2)*Power(y,2) + 1856*Power(y,4)) - 
          1536*Power(M,6)*(285*Power(x,6) + 624*Power(x,4)*Power(y,2) + 
             401*Power(x,2)*Power(y,4) + 62*Power(y,6)) + 
          1024*Power(M,5)*Sqrt(Power(x,2) + Power(y,2))*
           (784*Power(x,6) + 2151*Power(x,4)*Power(y,2) + 
             1968*Power(x,2)*Power(y,4) + 601*Power(y,6))) + 
       2*Power(a,8)*(-567*Power(M,10) - 
          9240*Power(M,9)*Sqrt(Power(x,2) + Power(y,2)) - 
          3360*Power(M,8)*(19*Power(x,2) + 18*Power(y,2)) - 
          320*Power(M,7)*Sqrt(Power(x,2) + Power(y,2))*
           (751*Power(x,2) + 623*Power(y,2)) - 
          12288*Power(Power(x,2) + Power(y,2),3)*(-Power(x,4) + 9*Power(y,4)) + 
          2048*M*Power(Power(x,2) + Power(y,2),2.5)*
           (55*Power(x,4) + 66*Power(x,2)*Power(y,2) + 35*Power(y,4)) - 
          1536*Power(M,5)*Sqrt(Power(x,2) + Power(y,2))*
           (349*Power(x,4) + 435*Power(x,2)*Power(y,2) + 98*Power(y,4)) - 
          256*Power(M,2)*Power(Power(x,2) + Power(y,2),2)*
           (-1609*Power(x,4) - 626*Power(x,2)*Power(y,2) + 263*Power(y,4)) - 
          480*Power(M,6)*(1067*Power(x,4) + 1734*Power(x,2)*Power(y,2) + 
             675*Power(y,4)) + 1024*Power(M,3)*Sqrt(Power(x,2) + Power(y,2))*
           (487*Power(x,6) + 1089*Power(x,4)*Power(y,2) + 
             753*Power(x,2)*Power(y,4) + 151*Power(y,6)) + 
          256*Power(M,4)*(-95*Power(x,6) + 775*Power(x,4)*Power(y,2) + 
             1691*Power(x,2)*Power(y,4) + 821*Power(y,6))) + 
       4*Power(a,4)*(-81*Power(M,14) - 
          1848*Power(M,13)*Sqrt(Power(x,2) + Power(y,2)) + 
          98304*(-Power(x,2) + Power(y,2))*Power(Power(x,2) + Power(y,2),6) + 
          32768*M*Power(Power(x,2) + Power(y,2),5.5)*
           (-19*Power(x,2) + 2*Power(y,2)) - 
          336*Power(M,12)*(56*Power(x,2) + 55*Power(y,2)) - 
          96*Power(M,11)*Sqrt(Power(x,2) + Power(y,2))*
           (1163*Power(x,2) + 1099*Power(y,2)) - 
          4096*Power(M,2)*Power(Power(x,2) + Power(y,2),4)*
           (169*Power(x,4) + 186*Power(x,2)*Power(y,2) + 29*Power(y,4)) + 
          8192*Power(M,3)*Power(Power(x,2) + Power(y,2),3.5)*
           (251*Power(x,4) + 370*Power(x,2)*Power(y,2) + 107*Power(y,4)) + 
          4096*Power(M,4)*Power(Power(x,2) + Power(y,2),3)*
           (1615*Power(x,4) + 2453*Power(x,2)*Power(y,2) + 826*Power(y,4)) + 
          2048*Power(M,5)*Power(Power(x,2) + Power(y,2),2.5)*
           (3919*Power(x,4) + 6252*Power(x,2)*Power(y,2) + 2345*Power(y,4)) - 
          128*Power(M,9)*Sqrt(Power(x,2) + Power(y,2))*
           (7549*Power(x,4) + 13685*Power(x,2)*Power(y,2) + 6148*Power(y,4)) - 
          48*Power(M,10)*(8701*Power(x,4) + 16434*Power(x,2)*Power(y,2) + 
             7737*Power(y,4)) + 256*Power(M,6)*Power(Power(x,2) + Power(y,2),2)*
           (18673*Power(x,4) + 31538*Power(x,2)*Power(y,2) + 12985*Power(y,4)) + 
          1024*Power(M,7)*Sqrt(Power(x,2) + Power(y,2))*
           (679*Power(x,6) + 2019*Power(x,4)*Power(y,2) + 
             2007*Power(x,2)*Power(y,4) + 667*Power(y,6)) - 
          128*Power(M,8)*(8555*Power(x,6) + 23039*Power(x,4)*Power(y,2) + 
             20437*Power(x,2)*Power(y,4) + 5953*Power(y,6)))))/
   Power(Power(a,8) + Power(M,8) + 16*Power(M,7)*Sqrt(Power(x,2) + Power(y,2)) + 
     112*Power(M,6)*(Power(x,2) + Power(y,2)) + 
     448*Power(M,5)*Power(Power(x,2) + Power(y,2),1.5) + 
     1120*Power(M,4)*Power(Power(x,2) + Power(y,2),2) + 
     1792*Power(M,3)*Power(Power(x,2) + Power(y,2),2.5) + 
     1792*Power(M,2)*Power(Power(x,2) + Power(y,2),3) + 
     1024*M*Power(Power(x,2) + Power(y,2),3.5) + 
     256*Power(Power(x,2) + Power(y,2),4) - 
     4*Power(a,6)*(Power(M,2) - 4*Power(y,2) + 4*M*Sqrt(Power(x,2) + Power(y,2))) + 
     Power(a,4)*(6*Power(M,4) + 48*Power(M,3)*Sqrt(Power(x,2) + Power(y,2)) - 
        64*M*Power(Power(x,2) + Power(y,2),1.5) + 
        16*Power(M,2)*(7*Power(x,2) + 5*Power(y,2)) + 
        32*(-Power(x,4) + 2*Power(x,2)*Power(y,2) + 3*Power(y,4))) + 
     4*Power(a,2)*(-Power(M,6) - 12*Power(M,5)*Sqrt(Power(x,2) + Power(y,2)) - 
        96*Power(M,3)*Power(Power(x,2) + Power(y,2),1.5) + 
        64*Power(y,2)*Power(Power(x,2) + Power(y,2),2) + 
        64*M*Power(Power(x,2) + Power(y,2),2.5) - 
        4*Power(M,4)*(14*Power(x,2) + 13*Power(y,2)) - 
        16*Power(M,2)*(Power(x,4) + 4*Power(x,2)*Power(y,2) + 3*Power(y,4))),3);
  if (PetscAbsReal(value) < EPSILON) return (0.);
  else return (value);
}
static PetscReal barrier(PetscReal x, PetscReal y, 
                                PetscReal rhori) {
  PetscReal barrier, a, r1,r2;
  
  a = PetscPowReal(10.0,20.0);
  r1 = PetscSqrtReal((x - 0.001)*(x - 0.001) + (y - 3)*(y - 3)) - rhori;
  r2 = PetscSqrtReal((x - 0.001)*(x - 0.001) + (y + 3)*(y + 3)) - rhori;
  barrier = -1.5/PETSC_PI*PetscAtanReal(r1*a)
            - 1.5/PETSC_PI*PetscAtanReal(r2*a) + 1.5/2 + 1.5/2;
  if(PetscAbsReal(barrier) < 0.1) return(0.);
  else return(barrier);
}

static PetscReal barrier_rup(PetscReal x, PetscReal y, 
			   PetscReal r0, PetscReal epsilon) {
  PetscReal barrier_rup, a, r1, r2;

  a = PetscPowReal(10.0, 20.0);
  r1 = PetscSqrtReal((x - 0.001)*(x - 0.001) + (y - 3)*(y - 3)) 
		- (r0 - epsilon);
  r2 = PetscSqrtReal((x - 0.001)*(x - 0.001) + (y - 3)*(y - 3)) 
		- (r0 + epsilon);
  barrier_rup = 1.0/PETSC_PI*(PetscAtanReal(a*r1)
		- PetscAtanReal(a*r2)); 
  
  return(barrier_rup);
}

static PetscReal barrier_rdown(PetscReal x, PetscReal y,
                             PetscReal r0, PetscReal epsilon) {
  PetscReal barrier_rdown, a, r1, r2;

  a = PetscPowReal(10.0, 20.0);
  r1 = PetscSqrtReal((x - 0.001)*(x - 0.001) + (y + 3)*(y + 3)) 
		- (r0 - epsilon);
  r2 = PetscSqrtReal((x - 0.001)*(x - 0.001) + (y + 3)*(y + 3)) 
		- (r0 + epsilon);
  barrier_rdown = 1.0/PETSC_PI*(PetscAtanReal(a*r1)
                  - PetscAtanReal(a*r2));
  
  return(barrier_rdown);
} 

static PetscReal barrier_thetaup(PetscReal x, PetscReal y,
			       PetscReal theta0, PetscReal epsilon) {
  PetscReal barrier_thetaup, a, theta1, theta2;
  
  a = PetscPowReal(10.0, 20.0);
  theta1 = PetscAtanReal((y - 3.0)/x) - (theta0 - epsilon);
  theta2 = PetscAtanReal((y - 3.0)/x) - (theta0 + epsilon);
  barrier_thetaup = 1.0/PETSC_PI*(PetscAtanReal(a*theta1)
                - PetscAtanReal(a*theta2));
  
  return (barrier_thetaup);
}

static PetscReal barrier_thetadown(PetscReal x, PetscReal y,
                               PetscReal theta0, PetscReal epsilon) {
  PetscReal barrier_thetadown, a, theta1, theta2;

  a = PetscPowReal(10.0, 20.0);
  theta1 = PetscAtanReal((y + 3.0)/x) - (theta0 - epsilon);
  theta2 = PetscAtanReal((y + 3.0)/x) - (theta0 + epsilon);
  barrier_thetadown = 1.0/PETSC_PI*(PetscAtanReal(a*theta1)
                - PetscAtanReal(a*theta2));
  
  return (barrier_thetadown);
} 

static PetscReal qro2z2(PetscReal ro, PetscReal z,
                        PetscReal a, PetscReal sig1,
                        PetscReal sig2) {
  PetscReal exp1, part1, part2, part3, part4, qroz;

  exp1 = PetscExpReal(-ro*ro/(sig1*sig1) - z*z/(sig2*sig2));
  part1 = -6.0*a*ro*exp1/(sig1*sig1);
  part2 = -2.0*a*ro*exp1/(sig2*sig2);
  part3 = 4.0*a*ro*z*z*exp1/(sig2*sig2*sig2*sig2);
  part4 = 4.0*a*ro*ro*ro*exp1/(sig1*sig1*sig1*sig1);

  qroz = part1 + part2 + part3 + part4;
  if(PetscAbsReal(qroz) < EPSILON) return(0.);
  else return(qroz);
}

static PetscReal qro(PetscReal ro, PetscReal z,
		     PetscReal a, PetscReal sig1,
		     PetscReal sig2) {
  PetscReal qro, exp1;
  exp1 = PetscExpReal(-ro*ro/(sig1*sig1) - z*z/(sig2*sig2));
  qro = a*exp1 - 2.0*a*ro*ro*exp1/(sig1*sig1);

  return (qro);
}

static PetscReal qz(PetscReal ro, PetscReal z,
		    PetscReal a, PetscReal sig1,
		    PetscReal sig2) {
  PetscReal qz, exp1;
  exp1 = PetscExpReal(-ro*ro/(sig1*sig1) - z*z/(sig2*sig2));
  qz = -2.0*a*ro*z*exp1/(sig2*sig2);

  return (qz);
}

static PetscReal q_func(PetscReal ro, PetscReal z,
		 	PetscReal a, PetscReal sig1,
			PetscReal sig2) {
  PetscReal exp1, qq;
  exp1 = PetscExpReal(-ro*ro/(sig1*sig1) - z*z/(sig2*sig2));
  qq = a*ro*exp1;
  
  return (qq);
}
 
static PetscReal NN(PetscReal ro, PetscReal z) {
  PetscReal NN;

  NN = 1.0;
  return NN;
}

static PetscReal dNdro(PetscReal ro, PetscReal z) {
  PetscReal dNdro;

  dNdro = 0.0;
  return dNdro;
}

static PetscReal dNdz(PetscReal ro, PetscReal z) {
  PetscReal dNdz;

  dNdz = 0.0;
  return dNdz;
}
 
extern PetscErrorCode InitialState(DM, Vec,  AppCtx*);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo*, Field**,
                                           Field**,  AppCtx*);
extern PetscErrorCode indexfunc_up(DMDALocalInfo*, PetscReal, PetscReal,PetscReal,
				PetscReal*, Field**,  AppCtx*); 
extern PetscErrorCode indexfunc_down(DMDALocalInfo*, PetscReal, PetscReal,PetscReal,
                                PetscReal*, Field**,  AppCtx*);
extern PetscErrorCode AngMom_ADMMass(DMDALocalInfo*, DM, Vec, PetscReal, AppCtx*);


//STARTMAIN
int main(int argc,char **argv) {
  PetscErrorCode ierr;
  DM            da;
  SNES          snes;
  AppCtx        user;
  Vec           x, Yexact, err, err1, y;
  PetscReal     AngMom = 0.0, hx,hy, errnorm, errnorm1, errnorm2;
  DMDALocalInfo info;

  ierr = PetscInitialize(&argc,&argv,NULL,help); if (ierr) return ierr;
  user.a  = 4.0;
  user.M  = 1.0;
  user.L  = 0.001;
  user.R  = 40.0;
  user.U  = 10.0;
  user.D  = -10.0;
  user.sig1 = 1.0;
  user.sig2 = 1.0;
  user.coer = 0.5;
  user.coethe = 0.5;
 
  ierr = DMDACreate2d(PETSC_COMM_WORLD,
               DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
               DMDA_STENCIL_STAR,  // for 5-point stencil
               5, 5, PETSC_DECIDE,PETSC_DECIDE,
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

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetDM(snes,da); CHKERRQ(ierr);
  ierr = DMDASNESSetFunctionLocal(da,INSERT_VALUES,
             (DMDASNESFunction)FormFunctionLocal,&user); CHKERRQ(ierr);
  //ierr = DMDASNESSetJacobianLocal(da,
    //         (DMDASNESJacobian)FormJacobianLocal,&user); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
  
  ierr = SNESSolve(snes,NULL,x); CHKERRQ(ierr);
  ierr = SNESComputeFunction(snes, x, y); CHKERRQ(ierr);
  ierr = AngMom_ADMMass(&info, da, x, AngMom, &user);
 
  hx = (user.R - user.L) / (PetscReal)(info.mx - 1);
  hy = (user.U - user.D) / (PetscReal)(info.my - 1);
  ierr = VecNorm(y, NORM_INFINITY, &errnorm); CHKERRQ(ierr);
  ierr = VecNorm(y, NORM_1, &errnorm1); CHKERRQ(ierr);

  ierr = VecNorm(y, NORM_INFINITY, &errnorm2); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,
           "L-Inf(res)=%g, L1(res)=%g, L-Inf(res)=%g,Angmom = \n", errnorm,errnorm1*hx*hy,errnorm2, AngMom); CHKERRQ(ierr);

  VecDestroy(&x);    VecDestroy(&Yexact); 
  VecDestroy(&err);   VecDestroy(&err1);
  SNESDestroy(&snes);    DMDestroy(&da);
  return PetscFinalize();
}
//ENDMAIN

//Initial state
PetscErrorCode InitialState(DM da, Vec Y, AppCtx* user) {
  PetscErrorCode   ierr;
  DMDALocalInfo    info;
  PetscInt         i,j;
  PetscReal        x, y, hx, hy, r0;
  DMDACoor2d       **aC;
  Field        **aY;

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
          if (i==info.mx-1 || j==0 || j==info.my-1){
                r0 = PetscSqrtReal(x*x + y*y);
                aY[j][i].u = 0.0;
		        aY[j][i].v = 1.0/(r0*r0*r0);
             }
             else{
                aY[j][i].u = 1.0 + 0.0*(1.5 - barrier(x, y, 1.5 - 0.52*hx))
                                 + 0.8*barrier(x, y, 1.5 - 0.52*hx);
		        aY[j][i].v =  0.02*barrier(x, y, 1.5 - 0.52*hx);
             }
      }
  }
  ierr = DMDAVecRestoreArray(da,Y,&aY); CHKERRQ(ierr);
  ierr = DMDARestoreCoordinateArray(da,&aC); CHKERRQ(ierr);
  return 0;
}

//index function for y > 0 part.
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
        value[2] = aY[j][i].v;
		k = k + 1;
          }
	  else if (barrier2 > 0.9) {
		value[1] = aY[j][i].u;
        value[3] = aY[j][i].v;
		l = l + 1;
	  }
	  else if (k == 1 && l == 1) {
	    goto here;
	  }
      }
  }
  here:
  if (l == 0 || k == 0){
    PetscPrintf(PETSC_COMM_WORLD, "l = %g,k = %g\n", l,k);
  }
  return 0.;
}

//index function for y < 0 part.
PetscErrorCode indexfunc_down(DMDALocalInfo *info, PetscReal theta,
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
          barrier1 = barrier_rdown(x, y, r1, coer*hx)
                        *barrier_thetadown(x, y, theta, coethe*hx);
	  barrier2 = barrier_rdown(x, y, r2, coer*hx)
                         *barrier_thetadown(x, y, theta, coethe*hx);
          if (barrier1 > 0.9) {
                value[0] = aY[j][i].u;
                value[2] = aY[j][i].v;
		k = k + 1;
          }     
	  else if (barrier2 > 0.9) {
		value[1] = aY[j][i].u;
        value[3] = aY[j][i].v;
		l = l + 1;
	  }
	  else if (k == 1 && l == 1) {
		goto here;
	  }
      }   
  }   
  here:
  if (l == 0 || k == 0){
    PetscPrintf(PETSC_COMM_WORLD, "l = %g,k = %g\n", l,k);
  }
  return 0.;
}

PetscErrorCode FormFunctionLocal(DMDALocalInfo *info, Field **aY, Field **aG, AppCtx *user) {
  PetscErrorCode   ierr;
  PetscInt   i, j, mx = info->mx, my = info->my, n, k;
  PetscReal  x, y, barrier_rup1, value[4], uxx, uyy,
 	     uxl, vxl, vxx, vyy, u1, u2, u3, v1, v2, v3,
	     uxr, uyl, uyr, vxr, vyl, vyr, part, theta = 0.0,
         coer = user->coer, dOmega, AngMom = 0.0, psi_BL;
  PetscReal  hx = (user->R - user->L) / (PetscReal)(mx - 1),
	     hy = (user->U - user->D) / (PetscReal)(my - 1),
         a1 = 0.99, a2 =  0.99;

  for (j = info->ys; j < info->ys + info->ym; j++) {
      y = hy * j + user->D;
      for (i = info->xs; i < info->xs + info->xm; i++) {
	      x = hx * i + user->L;
          if(i == mx - 1 || j == 0 || j == my - 1) {
	          aG[j][i].u = 0.0;
	          aG[j][i].v = 0.0;
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
              
              psi_BL = 1.0 + PetscSqrtReal(1.0 - a1*a1)/(2.0*PetscSqrtReal(x*x + (y - 3)*(y - 3))) 
                           + PetscSqrtReal(1.0 - a2*a2)/(2.0*PetscSqrtReal(x*x + (y + 3)*(y + 3)));

              uxx = (uxl - 2.0 * aY[j][i].u + uxr)/(hx*hx);
              uyy = (uyl - 2.0 * aY[j][i].u + uyr)/(hy*hy);
              //u1 = (aY[j][i+1].u - uxl) / (2.0*hx*x);
              //u2 = 1.0/4.0*qro2z2(x,y,a,sig1,sig2)*aY[j][i].u;
              u3 = 1.0/16.0*x*x
                     * (PetscPowReal((dOmegadrho(x, y - 3.0, 1.0, a1) 
                                        + dOmegadrho(x, y + 3.0, 1.0, a2)
                                        + 1.0/(2.0*hx)*(aY[j][i+1].v - vxl)),2.0)
                       +PetscPowReal((dOmegadz(x, y - 3.0, 1.0, a1)
                                        + dOmegadz(x, y + 3.0, 1.0, a2)
                                        + 1.0/(2.0*hy)*(aY[j+1][i].v - aY[j-1][i].v)), 2.0)) 
                     * PetscPowReal(aY[j][i].u + psi_BL, -7.0);
                //        + PetscPowReal(dOmegadz(x, y - 3.0, 1.0, a0),2.0)
                //        + PetscPowReal(dOmegadrho(x, y + 3.0, 1.0, -a0),2.0)
                //        + PetscPowReal(dOmegadz(x, y + 3.0, 1.0, -a0),2.0)
                //        + PetscPowReal(1.0/(2.0*hx)*(aY[j][i+1].v - vxl),2.0)
                //        + PetscPowReal(1.0/(2.0*hy)*(aY[j+1][i].v - aY[j-1][i].v),2.0))
                //     * PetscPowReal(aY[j][i].u + psi_BL, -7.0);

              vxx = (vxl - 2.0 * aY[j][i].v + vxr)/(hx*hx);
              vyy = (vyl - 2.0 * aY[j][i].v + vyr)/(hy*hy);
              v1 = 3.0/x*1.0/(2.0*hx)*(aY[j][i+1].v - vxl);
              //v2 = dNdro(x,y)*1.0/(2.0*hx)*(aY[j][i+1].v - vxl)
              //          + dNdz(x,y)*1.0/(2.0*hy)*(aY[j+1][i].v - aY[j-1][i].v);
              //v3 = 6.0/(aY[j][i].u)*1.0/(2.0*hx)*(aY[j][i+1].u - uxl)
              //                  *1.0/(2.0*hx)*(aY[j][i+1].v - vxl)
              //     + 6.0/(aY[j][i].u)*1.0/(2.0*hy)*(aY[j+1][i].u - aY[j-1][i].u)
              //                  *1.0/(2.0*hy)*(aY[j+1][i].v - aY[j-1][i].v);
	        
              v2 = LapOmega12(x, y - 3.0, 1.0, a1) + LapOmega12(x, y + 3.0, 1.0, a2);
              aG[j][i].u = 1.0/1.5*(uxx + uyy + u3);
                               //             *(1.5 - barrier(x,y,1.5-0.52*hx))
                               // + 0.0*barrier(x, y, 1.5-0.52*hx);
              aG[j][i].v = 1.0/1.5*(vxx + vyy + v1 - v2);
            }
      }
  }
  return 0.; 

}

PetscErrorCode AngMom_ADMMass(DMDALocalInfo *info, DM da, Vec Y, 
                                PetscReal AngMom, AppCtx *user) {
  PetscErrorCode   ierr;
  PetscInt         i, j, n, k = 0;
  PetscReal        x, y, hx, hy, dOmega, barrier_rup1,
                   barrier_thetaup1,theta = 0.0,Mass=0.0,
                   coer = user->coer, value[4],
                   mx = info->mx, my = info->my,
                   a = user->a, sig1 = user->sig1,
                   sig2 = user->sig2,barrier_thetadown1;
  DMDACoor2d       **aC;
  Field            **aY;

  ierr = DMDAGetCoordinateArray(da,&aC); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,Y,&aY); CHKERRQ(ierr);
  hx = (user->R - user->L) / (PetscReal)(info->mx - 1);
  hy = (user->U - user->D) / (PetscReal)(info->my - 1);
 
  n = (PetscInt)(PETSC_PI/(2.0*coer*hx));
  while (k < n) {
      theta = (k + 1)*2.0*coer*hx - PETSC_PI/2.0;
      for (j = info->ys; j < info->ys + info->ym; j++) {
          y = hy * j + user->D;
          for (i = info->xs; i < info->xs + info->xm; i++) {
	          x = hx * i + user->L;
	          barrier_rup1 = barrier_rdown(x, y, 1.5, coer*hx)
                                *barrier_thetadown(x, y, theta, coer*hx);
              if (barrier_rup1 > 0.6) {
                  ierr = indexfunc_down(info, theta, 1.5 + 2.0*coer*hx,
				        1.5 - 2.0*coer*hx, value, aY, user); CHKERRQ(ierr);
                  dOmega = (value[2] - aY[j][i].v)/(2.0*coer*hx);
                  //AngMom = AngMom + PetscAbsReal(PetscSinReal(theta))
                  //          *2.0*coer*hx; 
                  AngMom = AngMom + 1.0/8.0*1.5*1.5*2.0*coer*hx
                              * PetscPowReal(PetscAbsReal(PetscSinReal(theta)),3.0)
                              * PetscExpReal(q_func(x,y,a,sig1,sig2))*1.5*1.5
                              * dOmega*aY[j][i].u*aY[j][i].u;
                  //            * aY[j][i].u*aY[j][i].u*x*x*dOmega;
                  PetscPrintf(PETSC_COMM_WORLD,
                "AngMom=%g, dOmega=%g,theta=%g\n", AngMom, dOmega, theta); CHKERRQ(ierr);
                goto here;
               }
           }    
       }
       here:
       k = k + 1;
  }      
  for (j = info->ys; j < info->ys + info->ym; j++) {
          y = hy * j + user->D;
          for (i = info->xs; i < info->xs + info->xm; i++) {
	          x = hx * i + user->L;
             if (i == mx - 1){
                   Mass = Mass + (aY[j][i].u - aY[j][i-1].u)/hx*hy;
               }
               if (j == 0){
                   Mass = Mass + (aY[j][i].u - aY[j+1][i].u)/hy*hx;
               }
               if (j == my - 1){
                   Mass = Mass + (aY[j][i].u - aY[j-1][i].u)/hy*hx;
               }
 
          }
  }
  PetscPrintf(PETSC_COMM_WORLD,
           "AngMom=%g,Mass=%g\n", AngMom, Mass); CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(da,Y,&aY); CHKERRQ(ierr);
  ierr = DMDARestoreCoordinateArray(da,&aC); CHKERRQ(ierr);
  return 0;
}


