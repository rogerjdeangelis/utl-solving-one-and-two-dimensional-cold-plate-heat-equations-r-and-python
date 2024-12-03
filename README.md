# utl-solving-one-and-two-dimensional-cold-plate-heat-equations-r-and-python
Solving one and two dimensional cold plate heat equations r and python 
    %let pgm=utl-solving-one-and-two-dimensional-cold-plate-heat-equations-r-and-python;

    %stop_submission;

    I am a little outside my comfort zone here.

    Solving one and two dimensional cold plate heat equations r and python

    Graphic output 2 dimensional plate
    https://tinyurl.com/5n97pucx
    https://github.com/rogerjdeangelis/utl-solving-one-and-two-dimensional-cold-plate-heat-equations-r-and-python/blob/main/heatequ.png

    ASCII PLOT

                   X                  Symbol
         0        10        20
    Y ---+---------+---------+--- Y   .....   0.0 - 0.2 ring 1
      |                         |     +++++   0.2 - 0.4 ring 2
      |                         |     .....   0.4 - 0.6 ring 3
    20+       ....++++....      + 20  +++++   0.6 - 0.8 ring 4
      |    ..++++......++++..   |     .....   0.8 - 1.0 ring 5
      |   ..++..++++++++..++..  |
      |   ..+..++......++..+..  |
      |   ..+.++........++.+..  |
    10+   ..+..++......++..+..  + 10
      |   ..++.++++..++++.++..  |
      |    ..++....++....++..   |
      |      ..++.....+++..     |
      |       ....++++....     .|
     0+                         +  0
      |                         |
      ---+---------+---------+---
         0        10        20
                  X

    github
    https://tinyurl.com/2dxceary
    https://github.com/rogerjdeangelis/utl-solving-one-and-two-dimensional-cold-plate-heat-equations-r-and-python

          SOLUTIONS

              1. Two dimensional cold plate temperature distribution

              2. Python 1 dimensional closed form solution.
                 Temp(x) along the length of the cold plate
                 Requires boundary temperatures
                 (for symbolic mathematics I like sympy)

                 Solution: (given boundary conditions)
                                                          P*x
                                               --------------------------
                                               L*c_p*m_dot*(T_in - T_out)
                 T(x)= T_out + (T_in - T_out)*e

              3. Other sympy repos

    /*          ____     _             _     _       _       _
    / |  _ __  |___ \ __| |   ___ ___ | | __| |_ __ | | __ _| |_ ___
    | | | `__|   __) / _` |  / __/ _ \| |/ _` | `_ \| |/ _` | __/ _ \
    | | | |     / __/ (_| | | (_| (_) | | (_| | |_) | | (_| | ||  __/
    |_| |_|    |_____\__,_|  \___\___/|_|\__,_| .__/|_|\__,_|\__\___|
                                              |_|
    */

       Explanation  (2d temperture distribution)
       =========================================

         1. Discretization: The spatial domain is discretized into a grid with
            Nx and Ny points along the x and y directions respectively.
         2. Initial Condition: We initialize the temperature distribution using a given function f(x,y)
            f(x,y)  = sin(pi * x / Lx) * sin(pi * y / Ly
         3. Time-stepping: The code iteratively updates the temperature distribution using an explicit finite difference scheme.
         4. Visualization: The final temperature distribution is plotted using an image plot.

       This approach provides a numerical solution to the heat equation for a cold plate under specified conditions.
       Adjust parameters such as grid size and time step as needed for your specific application.


       INPUT: Parameters
       Lx <- 1       # Length of the plate in x-direction
       Ly <- 1       # Length of the plate in y-direction
       Nx <- 20      # Number of grid points in x-direction
       Ny <- 20      # Number of grid points in y-direction
       dx <- Lx / (Nx - 1)
       dy <- Ly / (Ny - 1)
       alpha <- 0.01 # Thermal diffusivity
       dt <- 0.001   # Time step
       Nt <- 100     # Number of time steps

       INPUT: Initial temperaure distribution
       f(x,y)  = sin(pi * x / Lx) * sin(pi * y / Ly

    /*
     _ __   _ __  _ __ ___   __ _ _ __ __ _ _ __ ___
    | `__| | `_ \| `__/ _ \ / _` | `__/ _` | `_ ` _ \
    | |    | |_) | | | (_) | (_| | | | (_| | | | | | |
    |_|    | .__/|_|  \___/ \__, |_|  \__,_|_| |_| |_|
           |_|              |___/
    */

    %utl_rbeginx;
    parmcards4;
    # Load necessary library
    library(pracma)
    library(haven)
    source("c:/oto/fn_tosas9x.R")
    # Parameters
    Lx <- 1       # Length of the plate in x-direction
    Ly <- 1       # Length of the plate in y-direction
    Nx <- 20      # Number of grid points in x-direction
    Ny <- 20      # Number of grid points in y-direction
    dx <- Lx / (Nx - 1)
    dy <- Ly / (Ny - 1)
    alpha <- 0.01 # Thermal diffusivity
    dt <- 0.001   # Time step
    Nt <- 100     # Number of time steps

    # Initialize temperature grid
    u <- matrix(0, nrow = Nx, ncol = Ny)

    # Initial condition: some function f(x,y)
    f <- function(x, y) {
      return(sin(pi * x / Lx) * sin(pi * y / Ly))
    }

    # Apply initial condition
    for (i in seq_len(Nx)) {
      for (j in seq_len(Ny)) {
        x <- (i - 1) * dx
        y <- (j - 1) * dy
        u[i, j] <- f(x, y)
      }
    }

    # Time-stepping loop
    for (n in seq_len(Nt)) {
      u_new <- u
      for (i in 2:(Nx-1)) {
        for (j in 2:(Ny-1)) {
          u_new[i, j] <- u[i, j] + alpha * dt * (
            (u[i+1, j] - 2*u[i, j] + u[i-1, j]) / dx^2 +
            (u[i, j+1] - 2*u[i, j] + u[i, j-1]) / dy^2
          )
        }
      }
      u <- u_new
    }
    str(u)
    png("d:/png/heatequ.png")
    # Plot the final temperature distribution
    image(u, main="Temperature Distribution", xlab="X", ylab="Y", col=heat.colors(100))
    png()
    xy=as.data.frame(u)
    fn_tosas9x(
          inp    = xy
         ,outlib ="d:/sd1/"
         ,outdsn ="want"
         )
    ;;;;
    %utl_rendx;

    proc print data=sd1.want;
    format _numeric_ 4.2;
    run;quit

    proc transpose data=sd1.want out=wantxpo;
    by rownames;
    run;quit;

    data wantxpox;
     set wantxpo;
     y=input(substr(_name_,2),3.);
     x=rownames;
     z=col1;
     drop _name_ col1;
    run;quit;

    options ls=64 ps=32;
    proc plot data=wantxpox(rename=y=y12345678901234567890);
    plot y12345678901234567890*x=z / contour=5  box;
    run;quit;

    /*           _               _
      ___  _   _| |_ _ __  _   _| |_
     / _ \| | | | __| `_ \| | | | __|
    | (_) | |_| | |_| |_) | |_| | |_
     \___/ \__,_|\__| .__/ \__,_|\__|
                    |_|
    */;

    /**************************************************************************************************************************/
    /*                                                                                                                        */
    /*   GRAPHIC                                                                                                              */
    /*   =======                                                                                                              */
    /*                                                                                                                        */
    /*   d:/png/heatequ.png                                                                                                   */
    /*                                                                            X                  Symbol          Z        */
    /*   HEAT MATRIX                                                    0        10        20                                 */
    /*   ===========                                               Y ---+---------+---------+--- Y   .....   0.0 - 0.2 ring 1 */
    /*                                                               |                         |     +++++   0.2 - 0.4 ring 2 */
    /*      V1      V2      V3 ...   V17     V18     V19             |                         |     .....   0.4 - 0.6 ring 3 */
    /*                                                             20+       ....++++....      + 20  +++++   0.6 - 0.8 ring 4 */
    /*  1  0.00    0.00    0.00 ...  0.00    0.00    0.00            |    ..++++......++++..   |     .....   0.8 - 1.0 ring 5 */
    /*  2  0.00    0.03    0.05 ...  0.08    0.05    0.03            |   ..++..++++++++..++..  |                              */
    /*  3  0.00    0.05    0.10 ...  0.15    0.10    0.05            |   ..+..++......++..+..  |                              */
    /*  ...                                                          |   ..+.++........++.+..  |                              */
    /* 18  0.00    0.05    0.10 ...  0.10    0.05    0.00          10+   ..+..++......++..+..  + 10                           */
    /* 19  0.00    0.03    0.05 ...  0.05    0.03    0.00            |   ..++.++++..++++.++..  |                              */
    /* 20  0.00    0.00    0.00 ...  0.00    0.00    0.00            |    ..++....++....++..   |                              */
    /*                                                               |      ..++.....+++..     |                              */
    /*                                                               |       ....++++....     .|                              */
    /*                                                              0+                         +  0                           */
    /*                                                               |                         |                              */
    /*                                                               ---+---------+---------+---                              */
    /*                                                                  0        10        20                                 */
    /*                                                                           X                                            */
    /*                                                                                                                        */
    /*                                                                                                                        */
    /*                                                                                                                        */
    /**************************************************************************************************************************/

    /*               _   _                   _     _        _                    _   __
    / |  _ __  _   _| |_| |__   ___  _ __   / | __| |   ___| | ___  ___  ___  __| | / _| ___  _ __ _ __ ___
    | | | `_ \| | | | __| `_ \ / _ \| `_ \  | |/ _` |  / __| |/ _ \/ __|/ _ \/ _` || |_ / _ \| `__| `_ ` _ \
    | | | |_) | |_| | |_| | | | (_) | | | | | | (_| | | (__| | (_) \__ \  __/ (_| ||  _| (_) | |  | | | | | |
    |_| | .__/ \__, |\__|_| |_|\___/|_| |_| |_|\__,_|  \___|_|\___/|___/\___|\__,_||_|  \___/|_|  |_| |_| |_|
     _  |_|    |___/     _
    (_)_ __  _ __  _   _| |_
    | | `_ \| `_ \| | | | __|
    | | | | | |_) | |_| | |_
    |_|_| |_| .__/ \__,_|\__|
            |_|
    */

       EXPLANATION

       One dimensional heat distribution along a length of coldplate

               /    -T_in + T(x) \
             P*|1 - -------------|
       dT      \    -T_in + T_out/
       -- =  ---------------------
       dx         L*c_p*m_dot

       Initial Condition

       T_in
       T_out
       L
       m_dot
       c_p
       P

       EASY TO SOLVE BECAUSE
       ======================

          dt
          -- = t(x)   (everthing else in the heat equation are constants)
          dx
               1
        Then  --- dt = dx
               t

        Integrate both sides

           / 1       /
           \ - dt  = \ dx
           / t       /

          ln(t)    =  x + c (c=constent of integration)

          Lets exponentiate both sides

                          x
             t(x)  = c * e

       LETS LET SYMPY DO THE MESSY ALGEBRA
       ====================================

                                          P*x
                               --------------------------
                               L*c_p*m_dot*(T_in - T_out)
       T_out + (T_in - T_out)*e

    /*           _   _
     _ __  _   _| |_| |__   ___  _ __    _ __  _ __ ___   __ _ _ __ __ _ _ __ ___
    | `_ \| | | | __| `_ \ / _ \| `_ \  | `_ \| `__/ _ \ / _` | `__/ _` | `_ ` _ \
    | |_) | |_| | |_| | | | (_) | | | | | |_) | | | (_) | (_| | | | (_| | | | | | |
    | .__/ \__, |\__|_| |_|\___/|_| |_| | .__/|_|  \___/ \__, |_|  \__,_|_| |_| |_|
    |_|    |___/                        |_|              |___/
    */

    /*--- Solving the first order differential equation of the heat equation (over time)   ---*/
    /*--- //www.perplexity.ai/search/plese-provide-a-solution-to-th-oPb23APfTTavhUfJ0cbNSw ---*/

    %utl_pybeginx;
    parmcards4;
    import sympy as sp
    from sympy import symbols, exp, pi, sqrt, integrate, diff, simplify, pprint, erf

    # Define symbols
    x = sp.Symbol('x')
    T = sp.Function('T')(x)
    T_in, T_out, L, m_dot, c_p, P, h = sp.symbols('T_in T_out L m_dot c_p P h')

    # Define the differential equation
    dT_dx = (P / (m_dot * c_p * L)) * (1 - (T - T_in) / (T_out - T_in))
    sp,pprint(dT_dx)
    eq = sp.Eq(T.diff(x), dT_dx)

    # Solve the differential equation
    solution = sp.dsolve(eq, T, ics={T.subs(x, 0): T_in})

    # Simplify the solution
    simplified_solution = sp.simplify(solution.rhs)

    print("Solution:")
    sp.pprint(simplified_solution)
    ;;;;
    %utl_pyendx;

    /*           _               _
      ___  _   _| |_ _ __  _   _| |_
     / _ \| | | | __| `_ \| | | | __|
    | (_) | |_| | |_| |_) | |_| | |_
     \___/ \__,_|\__| .__/ \__,_|\__|
                    |_|
    */

    /**************************************************************************************************************************/
    /*                                                                                                                        */
    /*  Solution:                                                                                                             */
    /*                                             P*x                                                                        */
    /*                                  --------------------------                                                            */
    /*                                  L*c_p*m_dot*(T_in - T_out)                                                            */
    /*    T(x)= T_out + (T_in - T_out)*e                                                                                      */
    /*                                                                                                                        */
    /**************************************************************************************************************************/

    /*____         _   _
    |___ /    ___ | |_| |__   ___ _ __   ___ _   _ _ __ ___  _ __  _   _  _ __ ___ _ __   ___  ___
      |_ \   / _ \| __| `_ \ / _ \ `__| / __| | | | `_ ` _ \| `_ \| | | || `__/ _ \ `_ \ / _ \/ __|
     ___) | | (_) | |_| | | |  __/ |    \__ \ |_| | | | | | | |_) | |_| || | |  __/ |_) | (_) \__ \
    |____/   \___/ \__|_| |_|\___|_|    |___/\__, |_| |_| |_| .__/ \__, ||_|  \___| .__/ \___/|___/
                                             |___/          |_|    |___/          |_|
    */

    REPO
    -------------------------------------------------------------------------------------------------------------------------------------
    https://github.com/rogerjdeangelis/utl-calculating-the-cube-root-of-minus-one-with-drop-down-to-python-symbolic-math-sympy
    https://github.com/rogerjdeangelis/utl-distance-between-a-point-and-curve-in-sql-and-wps-pythony-r-sympy
    https://github.com/rogerjdeangelis/utl-fun-with-sympy-infinite-series-and-integrals-to-define-common-functions-and-constants
    https://github.com/rogerjdeangelis/utl-maximum-likelihood-estimate-of--therate-parameter-lamda-of-a-Poisson-distribution-sympy
    https://github.com/rogerjdeangelis/utl-maximum-liklihood-regresssion-wps-python-sympy
    https://github.com/rogerjdeangelis/utl-mle-symbolic-solution-for-mu-and-sigma-of-normal-pdf-using-sympy
    https://github.com/rogerjdeangelis/utl-python-sympy-projection-of-the-intersection-of-two-parabolic-surfaces-onto-the-xy-plane-AI
    https://github.com/rogerjdeangelis/utl-r-python-compute-the-area-between-two-curves-AI-sympy-trapezoid
    https://github.com/rogerjdeangelis/utl-roots-of-a-non-linear-function-using-python-sympy
    https://github.com/rogerjdeangelis/utl-solve-a-system-of-simutaneous-equations-r-python-sympy
    https://github.com/rogerjdeangelis/utl-symbolic-algebraic-simplification-of-a-polynomial-expressions-sympy
    https://github.com/rogerjdeangelis/utl-symbolic-solution-for-the-gradient-of-the-cumulative-bivariate-normal-using-erf-and-sympy
    https://github.com/rogerjdeangelis/utl-symbolically-solve-for-the-mean-and-variance-of-normal-density-using-expected-values-in-SymPy
    https://github.com/rogerjdeangelis/utl-sympy-exact-pdf-and-cdf-for-the-correlation-coefficient-given-bivariate-normals
    https://github.com/rogerjdeangelis/utl-sympy-technique-for-symbolic-integration-of-bivariate-density-function
    https://github.com/rogerjdeangelis/utl-vertical-distance-covered-by-a-bouncing-ball-for-infinite-number-of-bounces-using-sympy

    https://github.com/rogerjdeangelis/utl-heat-map-of-correlation-matrix
    https://github.com/rogerjdeangelis/utl-three-heat-maps-of-bivariate-normal-wps-r-graph-plot
    https://github.com/rogerjdeangelis/utl_heatmap

    /*              _
      ___ _ __   __| |
     / _ \ `_ \ / _` |
    |  __/ | | | (_| |
     \___|_| |_|\__,_|

    */
