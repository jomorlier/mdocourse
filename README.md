# mdocourse 2020

****
LECTURES
****

The slides of the different lectures are available on *SUPAERO's LMS*. 

On some vulgarization articles

[Linkedin on mdo](https://www.linkedin.com/pulse/optimization-mdo-connecting-people-joseph-morlier/)

[Linkedin on mdo, but in french](https://www.linkedin.com/pulse/loptimisation-multidisciplinaire-pour-connecter-les-humains-morlier/)

[Theconversation on topopt, but also in french](http://theconversation.com/construire-une-aile-davion-en-lego-cest-possible-87126)


## Q&A session 
5 videoconferences (ZOOM) planned on:

1. Monday 20th of april 2020 9H30-10H30 Kick Off + help on assignement 1 fmincon(ZOOM, J. MORLIER) 
2. Wednesday 22th of april 2020 9H30-10H30 Q&A session and feedback on assignement 2 sensitivities (ZOOM, J. MORLIER) 
3. Monday 27th of april 2020 8H00-10H00  Q&A session and feedback on Exercice 1 truss optimization (ZOOM, M. HERRAZ, G. CAPASSO)
4. Monday 27th of april 2020 10H00-11H00  Q&A session and help on mini project topopt 3* (ZOOM, J. MORLIER) 
5. Wednesday 29th of april 2020 13H45-15H45  Q&A session and feedback on mini project topopt 3* (ZOOM, M. HERRAZ, G. CAPASSO) 


## AUDIO 

[Amanote's tutorial](https://www.youtube.com/watch?v=DvLyo9mtf3U)

[Course 1](https://github.com/jomorlier/mdocourse/blob/master/Course1.md)

[Course 1 outro](https://github.com/jomorlier/mdocourse/blob/master/Course1o.md)

[Course 2](https://github.com/jomorlier/mdocourse/blob/master/Course2.md)

[Course 3](https://github.com/jomorlier/mdocourse/blob/master/Course3.md)

****
TUTORIALS on GRADIENT/SENSITIVITY C1/C2
****


1st MATLAB tutorial on [Complex-step derivatives](http://htmlpreview.github.io/?https://github.com/jomorlier/mdocourse/blob/master/ComplexStep/ComplexStep.html)

Estimate derivatives by simply passing in a complex number to your function.
A single (complex) function evaluations computes both the function's value **(Re)** and the derivative **(Im)**.
Is it **always** possible to do this? I mean with a standard code form industry (Nastran, Fluent etc...)?

2nd MATLAB tutorial on [gradient evaluation](http://htmlpreview.github.io/?https://github.com/jomorlier/mdocourse/blob/master/Sensibility/sensitivity_TD.html)

Comparison of **Symbolic/Finite Differences/DIRECT/ADJOINT Method** on a really simple mechanical system (2DOFs).
Play with the code for **checking** Symbolic with Finite Differences. Play with $\Delta_x$ ?
By the way, just add the complex step approach, not so difficult when you have access to the **original** code.
Oh, at the end which method is exact? 

How Nastran is doing for gradient computation on SOL2OO
[gradient nastran](https://app.amanote.com/note-taking/document/827200fd-e137-475b-aab5-58d734086654)

****************

Tutorials on TOPOPT C3

****************





****************

ADVANCED TOPOPT (for interested MAE student)

****************


https://github.com/jomorlier/ALMcourse

Recent advances in Topology optimization (ALM) 23/11/18 [at ENSAM Bordeaux](https://github.com/mid2SUPAERO/Outputs/blob/master/Presentation_JMSC_FA.pdf)

You can find in top88 repository all the files, paper, and ...

[The 3 point bending projected corrected using top88](http://htmlpreview.github.io/?https://github.com/jomorlier/ALMcourse/blob/master/top88/topopt_3ptBENDING.html)

for people wondering about how linear elasticity of top88.m is working

Prof, can you help to us understand how to compute the stiffness matrix of 2D membrane?
Explictely ?  [Membrane2D_K](http://htmlpreview.github.io/?https://github.com/jomorlier/feacourse/blob/master/Membrane2D_K/Elementarystiffrecmesh.html)


Have a look in the top3D repository to do the same but in 3D !!! 


MATLAB tutorial Advanced Topology Optimization.

For people who are not familliar with Finite Element Method:
http://designinformaticslab.github.io/mechdesign_lecture/2018/01/28/featop.html

Part A:  [Constraints Agreggation](http://htmlpreview.github.io/?https://github.com/jomorlier/mdocourse/blob/master/AdvancedTopOpt/ConstraintsAgreggation.html)
Thanks to my PhD Simone.

Part B:  [Stress Based TopOpt](http://htmlpreview.github.io/?https://github.com/jomorlier/mdocourse/blob/master/AdvancedTopOpt/StressBasedTopOpt.html)
Thanks **AGAIN** to my PhD Simone.
Before you can use Method of Moving Asymptotes (MMA) as an optimizer in our stress based topology optimization program, you need to obtain the Matlab implementation of MMA from Prof. Krister Svanberg (krille@math.kth.se) from KTH in Stockholm Sweden.

****************

Gaussian Processes aka Kriging

****************

4th MATLAB tutorial on regression using [GP, or Kriging](http://htmlpreview.github.io/?https://github.com/jomorlier/mdocourse/blob/master/GP_Tutorial/GP_Tutorial.html)

It's always working ? On interpolation **and** extrapolation ?
By the way, GP is an interpolating **or** regressing method ? Play with the code.

My first book on  Kriging: [Aerospace Engineering](https://optimizationcodes.wordpress.com)
My first book on  GP:[Machine Learning](http://www.gaussianprocess.org/gpml/code/matlab/doc/)


Small recipes on Surrogate modeling and Bayesian optimization 3/10/19 [at TU DELFT](https://github.com/mid2SUPAERO/Outputs/blob/master/Recipes_DELFT-3-10-19-compressed.pdf)

Try our SMT Python package for creating Surrogate models

https://smt.readthedocs.io/en/latest/_src_docs/surrogate_models.html


Try our SMT Python package for doing Efficient Global Optimization 

https://smt.readthedocs.io/en/latest/_src_docs/applications/ego.html


nice application http://webfoil.engin.umich.edu/about

****************

ROM 

****************

Inspired from Menke's book: Another MATLAB tutorial on [Reduced Order Modeling](http://htmlpreview.github.io/?https://github.com/jomorlier/mdocourse/blob/master/RoM/ROM.html) using SVD/POD approach.

Sparse and Distributed Gaussian process for Flight test and Structural dynamics 12/06/17 [at 3AF BIG DATA](https://github.com/mid2SUPAERO/Outputs/blob/master/MDO_12-06-17_3AFBigData.pdf)

****************

MDA

****************
What about MDO ? No MDO without MDA...
Here is a simple example using 
[Fixed Point Iteration ](http://htmlpreview.github.io/?https://github.com/jomorlier/mdocourse/blob/master/MDA/tutorialFPI.html) 

****************

MDO

****************

A tiny introduction to MDO 21/2/19 [at ICA](https://github.com/mid2SUPAERO/Outputs/blob/master/Presentation_JM_MDO-compressed.pdf)

MDO in academia: From Classrooms to Research (and vice versa) 19/11/19 [at 2nd European Workshop on MDO for Industrial Applications in Aeronautics IRT SE](https://github.com/mid2SUPAERO/Outputs/blob/master/MDOinACADEMIA-compressed.pdf)

How do we use OpenMDAO in our Research activities at ONERA/SUPAERO (and also in classrooms)? 21/10/19 [at Ohio Aerospace Institute](https://github.com/mid2SUPAERO/Outputs/blob/master/OpenMDAO_Cleveland_LIGHT2-compressed.pdf)


## REFERENCES

Over the year, Thanks to my colleagues of UoM, ONERA, Supaero and Airbus.

Books
"Principles of Optimal Design  Modeling and Computation”, P.Y. Papalambros & D.J. Wilde, Cambridge University Press
“Elements of Structural Optimization”, R.T. Haftka & Z. Gurdal, Kluwer 
"Engineering design via surrogate modelling: a practical guide", Forrester, Alexander, Andras Sobester, and Andy Keane, , John Wiley & Sons, 2008.

Online courses
SMO course, Haftka, Florida Institue of Technology
Optimisation structurale et multidisciplinaire P. Duysinx, ULB
Multidisciplinary System Design Optimization, Prof. Olivier de Weck, Prof. Karen Willcox MIT Course Number: ESD.77 / 16.888
Design Optimization Tools : Course Material Pascal Etman, Piet Schreurs, Rob Bastiaans EIT Eindhoven
Engineering Optimization, Fred van Keulen, TU DELFT
Topopt.dk, Topology Optimization, Ole Sigmund , Mechanical Engineering Technical University of Denmark (DTU)
Optimisation de structures, Gregoire Allaire, X
Engineering problems in Matlab, Magrab
Optimisation numerique, Aleks Janka, Université de Fribourg



