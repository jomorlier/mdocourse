%% 
%% TP MDO - Part 2
%% Simone Coniglio - Joseph Morlier 
%% *Exercice 2 - Towards stress based topology optimization*
% In this exercice we will introduce stress based constraint in the classic 
% topology optimization framework of top88.m. You can also use top99.m code, the 
% only difference will be the computational efficiency in Matlab. 
% 
% In the first part of this exercice you will be guided to compute Von Mises 
% stress  in the same framework (top88.m). 
% 
% We will also treat the efficient implementation of adjoint approach for 
% constraint gradient evaluation of constraint gradient in a Finite Element framework.
% 
%  Finally we will use an the most used optimization solver in structural 
% optimization called MMA [3]. The following codes will be provided (mmasub.m, 
% subsolve.m, topstressMMA.m,). Why ? 
% 
% Because classical methods (interior point and sequential quadratic programming) 
% already implemented in fmincon Matlab will have a slow convergence rate due 
% to Hessian sparsity. 
%% Review of stress evaluation in 2D stress plane quadrilateral elements:
% In top88 and top99 we are using structured uniform mesh with constant 1x1 
% size, in the follow we get the stress for a general square element of side $2a\times2a$ 
% in the figure below: 
% 
% 
% 
% 
% 
% The shape functions $N_i(x,y)$ for $i=1,2,3,4$ are:
% 
% $$ N_i(x,y)=\frac{(1+\frac{x_i}{a^2}x)(1+\frac{y_i}{a^2}y)}{4} $$
% 
% So that the displacement field can be described as:
% 
% $$u(x,y)=\sum_ {i=1}^4N_i(x,y)u_i\\v(x,y)=\sum_ {i=1}^4N_i(x,y)v_i$$
% 
% The strain field is simply:
% 
%  $$\epsilon_x(x,y)=\frac{\partial u}{\partial x}(x,y)\\\epsilon_y(x,y)=\frac{\partial 
% v}{\partial y}(x,y)\\\gamma_{xy}(x,y)=\frac{\partial v}{\partial x}(x,y)+\frac{\partial 
% u}{\partial v}(x,y)$$
% 
%  For clarity we will compute the stress only in a gauss point located in 
% $(0,0)$. In this point the shape function derivatives are:
% 
% $$\frac{\partial N_i}{\partial x}(0,0)=\frac{x_i}{4a^2}\\\frac{\partial 
% N_i}{\partial y}(0,0)=\frac{y_i}{4a^2}$$
% 
% Then defining the dispacement vector $\lbrace q\rbrace$ as:
% 
% $$\lbrace q\rbrace=\left\lbrace u_1,v_1,u_2,v_2,u_3,v_3,u_4,v_4 \right\rbrace^T$$
% 
% with $u_1$etc... degrees of freedom depicted in the above picture.
% 
% The displacement to strain relationship can be defined as:
% 
%  $$\left\lbrace \begin{array}{c}\epsilon_x\\\epsilon_y\\\gamma_{xy}\end{array}\right\rbrace(0,0)=\frac{1}{4a}\left[\begin{array}{c} 
% -1&0&1&0&1&0&-1&0 \\ 0&-1&0&-1&0&1&0&1 \\-1&-1&-1&1&1&1&1&-1\end{array}\right]\lbrace 
% q\rbrace=[B]\lbrace q\rbrace$$
% 
% Once that the strain are linked to displacement the stress vector can also 
% be evauated by the Hooke's generalized law for planar stress defining the stress 
% to strain matrix $[D]$ as:
% 
% $$[D]=\frac{E}{1-\nu^2}\left[\begin{array}{ccc}1&\nu&0\\\nu&1&0\\0&0&\frac{1-\nu}{2}\end{array}\right]$$
% 
% so that 
% 
% $$\left\lbrace \begin{array}{c}\sigma_x\\\sigma_y\\\tau_{xy}\end{array}\right\rbrace(0,0)=[D]\left\lbrace 
% \begin{array}{c}\epsilon_x\\\epsilon_y\\\gamma_{xy}\end{array}\right\rbrace(0,0)=[D][B]\lbrace 
% q \rbrace$$
% 
% The Von Mises (VM) stress is the most common criteria used for metallic 
% material failure. Practicaly VM stress value can be seen as the pure normal 
% stress value that is as dangerous as the combined stress state given by the 
% stress tensor. In planar stress it is defined as:
% 
% $$\sigma_{VM}=\sqrt{\sigma_x^2+\sigma_y^2-\sigma_x\sigma_y+3\tau_{xy}^2}$$
% 
% that can be rewritten in matrix form as:
% 
% $$\sigma_{VM}=\sqrt{\left\lbrace\begin{array}{ccc}\sigma_x & \sigma_y&\tau_{xy}\end{array}\right\rbrace\left[\begin{array}{ccc}1& 
% -\frac{1}{2}&0\\-\frac{1}{2} &1 &0\\0&0&3\end{array}\right]\left\lbrace\begin{array}{c}\sigma_x 
% \\ \sigma_y \\ \tau_{xy}\end{array}\right\rbrace}=\sqrt{\lbrace q\rbrace^T [B]^T[D]\left[\begin{array}{ccc}1& 
% -\frac{1}{2}&0\\-\frac{1}{2} &1 &0\\0&0&3\end{array}\right][D][B]\lbrace q\rbrace}=\sqrt{\lbrace 
% q\rbrace^T [S]\lbrace q \rbrace}$$
% 
% Since our mesh is uniform $\frac{1}{E}\left[S\right]$ is the same for each 
% element and can be evalued once for all.
% 
% 
%% 2a) Verify the stress constraints in final solution of compliance based optimization (solution provided)
% Using the formulation above you should be able to evaluate the Von Mises stress 
% in each Element centroid. Modify slightly the 88 line Matlab code of Andreassen 
% including some line for the evaluation of $[S]$ and other for the stress. finally 
% use imagesc matlab function to get the a final stress plot of your design zone.What 
% do you remarck? Which is the final stress distribution? Can you say that your 
% desing is always feasible? 
% 
% *Solution (using compare button)*
% 
% In our code we made a few changes in the 88lines to fit the new specification 
% see the report of file comparison below:
% 
% 
% 
% With just these modifications we can get all the answers:

top88stress(100,50,0.1,3,2,2);
%% 
% 
% 
% Compliance driven topology optimization do not take into account for feasibility 
% constrants. The material could be failing when the volume fraction demanded 
% is too small.
% 
% The engineer can't know apriori which level of volume fraction will be 
% feasible. 
% 
% In the next we will implement a different approach that avoid this issue.
%% 2b) Solving stress based topology optimization problem using fmincon (on the MBB configuration)
% Given the design space $\mathbf{D}$ discretized with finite element we want 
% to find the best material layout described by the design variable vector $\lbrace 
% x \rbrace$, that minimize the employed mass  $M(\lbrace x \rbrace ) = \frac{\rho}{N}\sum_{i=1}^Nx_i$ 
% , still respecting stress constraints inside the design zone. For uniforme density 
% the mass is simply proportional to the volume fraction $V(\lbrace x \rbrace 
% ) = \frac{1}{N}\sum_{i=1}^Nx_i$ so that the stress based topology optimization 
% can be formulated as:
% 
% $$\left\lbrace\begin{array}{cc}\min_{\lbrace x\rbrace}V(\lbrace x\rbrace)\\g_i(\lbrace 
% u\left(\lbrace x \rbrace\right) \rbrace)=\frac{\sigma_i(\lbrace u\left(\lbrace 
% x \rbrace\right) \rbrace)}{\sigma_{lim}}-1\leq0 & \forall i\in\lbrace1,2,...,N\rbrace|x_i>0\\0\leq 
% x_i\leq 1 & \forall i\in\lbrace1,2,...,N\rbrace\\[K(\lbrace x \rbrace)]\lbrace 
% u\left(\lbrace x \rbrace\right) \rbrace=\lbrace F\rbrace\end{array}\right.$$
% 
% The last equation is the static balance of equation deriving from the Finite 
% Element Model of the design zone.
% 
% As for the two bar truss problem, this problem shows vanishing constraints. 
% The dumb solution $\lbrace x \rbrace=\lbrace 0 \rbrace$ is eliminated from the 
% design space due to the existance of a displacement solution.
% 
% The Young modulus is  obtained by the power low according to [].
% 
% $$E_i(x_i)=E_{min}+\left(E_{0}-E_{min} \right)x_i^p \quad p>1$$
% 
% To get access to singular optima that are also present in this problem, 
% the unified aggregation and relaxation technique of Verbart [1] is studied here.
% 
% Firstly the failure is associated not to the global stress (i.e. the V.M. 
% stress that we obtain using the real young module of the element) but to the 
% microscopic stress (i.e. the V.M. stress per density of material, see [6]). 
% The explaination given to this choice is that an element with an intermediate 
% density can be considered as a composite material layered with full material 
% and voids that have a thickness ratio of $\frac{x}{1-x}$. We get the microscopic 
% stress considering a $[D_0]$ matrix evaluated for full material module $E_0$.
% 
% We will note $\sigma_{i0}$ the microscopic stress derived by this definition 
% and as $g_{i0}$ the correspondinc constraint violation. We get the relaxed constraint 
% violation $\bar{g}_{i0}$simply multipling $g_{i0}$ by the corresponding density 
% $x_i$:
% 
% The new problem incorporating microscopic stress and relaxation is: 
% 
% $$\left\lbrace\begin{array}{cc}\min_{\lbrace x\rbrace}V(\lbrace x\rbrace)\\\bar{g}_i(\lbrace 
% x \rbrace,\lbrace u\left(\lbrace x \rbrace\right) \rbrace)=x_i\left(\frac{\sigma_{i0}(\lbrace 
% u\left(\lbrace x \rbrace\right) \rbrace)}{\sigma_{lim}}-1\right)\leq0 &\forall 
% i\in\lbrace1,2,...,N\rbrace\\0\leq x_i\leq 1 & \forall i\in\lbrace1,2,...,N\rbrace\\[K(\lbrace 
% x \rbrace)]\lbrace u\left(\lbrace x \rbrace\right) \rbrace=\lbrace F\rbrace\end{array}\right.$$
% 
% You could try to solve this problem with fmincon getting very slow in each 
% iteration due to the important number of non-linear constraints in your design 
% problem (please try :) ... ). To get faster in this case the aggregation is 
% needed. We will use again the lower bound KS-function to approximate the maximum 
% function.
% 
% $$\left\lbrace\begin{array}{cc}\min_{\lbrace x\rbrace}V(\lbrace x\rbrace)\\G_{KS}(\lbrace 
% x \rbrace,\lbrace u\left(\lbrace x \rbrace\right) \rbrace)=\frac{1}{P}\ln\left(\frac{1}{N}\sum_{i=1}^{N}e^{P\bar{g}_i(\lbrace 
% x \rbrace,\lbrace u\left(\lbrace x \rbrace\right) \rbrace)}\right)\leq0\\0\leq 
% x_i\leq 1 & \forall i\in\lbrace1,2,...,N\rbrace\\[K(\lbrace x \rbrace)]\lbrace 
% u\left(\lbrace x \rbrace\right) \rbrace=\lbrace F\rbrace\end{array}\right.$$
% 
% 
% 
% As in the two bar truss problem, you must create two functions: one for 
% the objective and one for the aggregated constraint.
% 
% The static balance equation will be embedded in the evaluation of $G_{KS}$ 
% function.
% 
% Try to solve the MMB problem with a very small mesh refinement (24*12). 
% Do you find an optimal solution? Can you detect a problem?
% 
% *Solution*
% 
% Finite difference approach is prohibitive in this case even for a coarse 
% mesh. Objective and constraint gradients must be provided.
%%  2c) The adjoint approach for stress gradient sensitivity
% The adjoint approach for the evaluation of the constraint is developped in 
% the follows in the general case of response function depending on FEA solution 
% displacement vector $\lbrace u\left(\lbrace x \rbrace\right) \rbrace$.
% 
% We want to evaluate the total response sensitivity the vector notation 
% $\lbrace\cdot\rbrace}$ will be neglected for brevity in this paragraph. 
% 
% $$\frac{d o(x,u(x))}{dx}=\frac{\partial o(x,u(x))}{\partial x}+\sum_{i=1}^{N_{DOF}}\frac{\partial 
% o(x,u(x))}{\partial u_j}\frac{\partial u_j(x)}{\partial x}$$
% 
% Displacement sensitivity are available using derving static balance equation:
% 
% $$\frac{d\left([K(\lbrace x \rbrace)]\lbrace u\left(\lbrace x \rbrace\right) 
% \rbrace\right)}{dx}=\frac{d\lbrace F\rbrace}{dx}=\lbrace 0\rbrace$$
% 
% So that 
% 
% $$\left[\frac{d \lbrace u\rbrace}{d \lbrace x\rbrace}\right]=-[K]^{-1}\left[\bar{\frac{dK}{dx}}\right]\lbrace 
% u \rbrace$$
% 
% Where the $\left[\bar{\frac{dK}{dx}}\right]$ is a third order array. The 
% above expression get easier considering one design variable at time.
% 
% $$\begin{array}{cccc}\frac{d o(x,u(x))}{dx_i}=\frac{\partial o(x,u(x))}{\partial 
% x_i}-&\left\lbrace\frac{\partial  o(x,u(x))}{\partial u}\right\rbrace^T&[K]^{-1}&\left[\frac{dK}{dx_i}\right]\lbrace 
% u\rbrace\\& 1\times N_{DOF}& N_{DOF}\times N_{DOF}&N_{DOF}\times 1\end{array}$$
% 
% From this equation it is clear that we have a choice, we can either compute 
% one inverse per each design variable computing the therm $[K]^{-1}\left[\frac{dK}{dx_i}\right]\lbrace 
% u\rbrace$ to get the final sensitivity, either compute just once for each response 
% the adjoint vector:
% 
% $$\lbrace \lambda\rbrace=[K]^{-1}\left\lbrace\frac{\partial  o(x,u(x))}{\partial 
% u}\right\rbrace$$
% 
% and than use this vector to compute all design variable sensitivities as:
% 
% $$\frac{d o(x,u(x))}{dx_i}=\frac{\partial o(x,u(x))}{\partial x_i}-\lbrace 
% \lambda\rbrace^T\left[\frac{dK}{dx_i}\right]\lbrace u\rbrace$$
% 
% Depending on the size of your problem (i.e. number of responses vs number 
% of variables) the first or the second strategy will be more efficient.
% 
% In the case of stress based topology optimization the number of design 
% variable is much greater then the number of response depending on displacement 
% field. Therefore the second approach will be preferred.
% 
% get the analytical expression of $\frac{\partial G_{KS}(x,u(x))}{\partial 
% x_i} $ and of $\left\lbrace\frac{\partial G_{KS}(x,u(x))}{\partial u} \right\rbrace$ 
% then compute the total sensitivity by adjoint approach. Provide this gradient 
% to the non-linear constraint function.
% 
% *Solution*
% 
% It can be showed that
% 
% $$\frac{\partial G_{KS}(x,u(x))}{\partial x_j} =\frac{\left(\frac{\sigma_{j0}}{\sigma_{lim}}-1\right)e^{P\bar{g}_j}}{\sum_{i=1}^{N}e^{P\bar{g}_i}}$$
% 
% $$\left\lbrace\frac{\partial G_{KS}(x,u(x))}{\partial u} \right\rbrace=\frac{\sum_{i=1}^{N}\frac{x_i}{\sigma_{lim}\sigma_{i0}}e^{P\bar{g}_i}[S_{0i}]\lbrace 
% u\rbrace}{\sum_{i=1}^{N}e^{P\bar{g}_i}}$$
% 
% **
% 
% Using the adjoint approach the gradient was provided to fmincon.
% 
% As a full hessian matrix is inversed at each optimization step, the computational 
% cost of this approach for topology optimization is still unfeasible!
% 
% try also the SQP algorithm. It is even worse!
% 
% 

topstress(24,12,3,2,2,1,4);
%% 
%  Note that the final desing is not even converged and has big grey area 
% that means either that the solver found a local minimum, either that the convergence 
% criteria could not be strong enough.
% 
% Looking at the microscopic stress in the element with a density greater 
% than 0.3 one can better understand why the optimization stopped. 
% 
% The stress is too big at the boundary condition, the solver stuck in a 
% local minimum.
% 
% One should try different mesh refinement in order to get a feasible. Using 
% fmincon this could be prohibitive.
% 
% That's why MMA is used in order to get more reasonable mesh refinement 
% in a finite amount of time.
% 
% Note that changing the mesh refinement and keeping the rest of the code 
% as it is the design zone change also its geometry. 
%%  2d (Bonus)) Try to use MMA optimization solver instead of fmincon
% You will be provided of MMA code and documentation of Svanberg. You can use 
% this code only citing his author Svanberg [3].
% 
% *solution*
% 
% This optimization solver is delicate and can be tuned refering to the reference 
% work of Svanberg and to the provided documentation.
% 
% MMA can carry out bigger optimization problem stil keeping fesible design.
% 
% The result is this time more understandable.
% 
% 

topstressMMA(50,25,3,2,2,1,4);
%% 
% 
% 
% Still we can note that the constraint at the BC is violated due to the 
% very small value of the aggregation parameter chosen in the unified approach.
% 
% On may also change the value of P in order to get feasible design.
% 
% Be aware that an increased value of P will also be the cause of a longer 
% optimization time.
% 
% 
%% References
% [1] Verbart, Alexander, Matthijs Langelaar, and Fred Van Keulen. "A unified 
% aggregation and relaxation approach for stress-constrained topology optimization." 
% _Structural and Multidisciplinary Optimization_ 55.2 (2017): 663-679.
% 
% [2] Coniglio, Simone, et al. "Original Pylon Architecture Design Using 
% 3D HPC Topology Optimization." _2018 AIAA/ASCE/AHS/ASC Structures, Structural 
% Dynamics, and Materials Conference_. 2018.
% 
% [3] Svanberg, K. (1987). The method of moving asymptotes—a new method for 
% structural optimization. _International journal for numerical methods in engineering_, 
% _24_(2), 359-373.
% 
% [4] Bendsøe, Martin P. "Optimal shape design as a material distribution 
% problem." _Structural optimization_ 1.4 (1989): 193-202.
% 
% [5] Zhou, M., and G. I. N. Rozvany. "The COC algorithm, Part II: topological, 
% geometrical and generalized shape optimization." _Computer Methods in Applied 
% Mechanics and Engineering_ 89.1-3 (1991): 309-336.
% 
% [6] Duysinx, Pierre, and Martin P. Bendsøe. "Topology optimization of continuum 
% structures with local stress constraints." _International journal for numerical 
% methods in engineering_ 43.8 (1998): 1453-1478.
% 
%