# mdo course 2020

Thanks to my colleagues of ONERA Toulouse

****************

LECTURES

****************


The slides of the different courses are available on *SUPAERO's LMS*. 

[Linkedin on mdo](https://www.linkedin.com/pulse/optimization-mdo-connecting-people-joseph-morlier/)

[Linkedin on mdo, but in french](https://www.linkedin.com/pulse/loptimisation-multidisciplinaire-pour-connecter-les-humains-morlier/)

[Theconversation on topopt, but also in french](http://theconversation.com/construire-une-aile-davion-en-lego-cest-possible-87126)


****************

TUTORIALS on GRADIENT/SENSITIVITY

****************


1st MATLAB tutorial on [Complex-step derivatives](http://htmlpreview.github.io/?https://github.com/jomorlier/mdocourse/blob/master/ComplexStep/ComplexStep.html)

Estimate derivatives by simply passing in a complex number to your function.
A single (complex) function evaluations computes both the function's value **(Re)** and the derivative **(Im)**.
Is it **always** possible to do this? I mean with a standard code form industry (Nastran, Fluent etc...)?

2nd MATLAB tutorial on [gradient evaluation](http://htmlpreview.github.io/?https://github.com/jomorlier/mdocourse/blob/master/Sensibility/sensitivity_TD.html)

Comparison of **Symbolic/Finite Differences/DIRECT/ADJOINT Method** on a really simple mechanical system (2DOFs).
Play with the code for **checking** Symbolic with Finite Differences. Play with $$\delta_x$$ ?
By the way, just add the complex step approach, not so difficult when you have access to the **original** code.
Oh, at the end which method is exact? 

How Nastran is doing for gradient computation on SOL2OO
[gradient nastran](https://app.amanote.com/note-taking/document/827200fd-e137-475b-aab5-58d734086654)

****************

TOPOPT

****************



https://github.com/jomorlier/ALMcourse

Recent advances in Topology optimization (ALM) 23/11/18 [at ENSAM Bordeaux](https://github.com/mid2SUPAERO/Outputs/blob/master/Presentation_JMSC_FA.pdf)

****************

ADVANCED TOPOPT

****************
3rd MATLAB tutorial Advanced Topology Optimization.

For people who are not familliar with Finite Element Method:
http://designinformaticslab.github.io/mechdesign_lecture/2018/01/28/featop.html

Part A:  [Constraints Agreggation](http://htmlpreview.github.io/?https://github.com/jomorlier/mdocourse/blob/master/AdvancedTopOpt/ConstraintsAgreggation.html)
Thanks to my PhD Simone.

Part B:  [Stress Based TopOpt](http://htmlpreview.github.io/?https://github.com/jomorlier/mdocourse/blob/master/AdvancedTopOpt/StressBasedTopOpt.html)
Thanks **AGAIN** to my PhD Simone.
Before you can use Method of Moving Asymptotes (MMA) as an optimizer in our stress based topology optimization program, you need to obtain the Matlab implementation of MMA from Prof. Krister Svanberg (krille@math.kth.se) from KTH in Stockholm Sweden.

****************

Gaussian Processes

****************

4th MATLAB tutorial on regression using [GP, or Kriging](http://htmlpreview.github.io/?https://github.com/jomorlier/mdocourse/blob/master/GP_Tutorial/GP_Tutorial.html)

It's always working ? On interpolation **and** extrapolation ?
By the way, GP is an interpolating **or** regressing method ? Play with the code.


Small recipes on Surrogate modeling and Bayesian optimization 3/10/19 [at TU DELFT](https://github.com/mid2SUPAERO/Outputs/blob/master/Recipes_DELFT-3-10-19-compressed.pdf)

Try our SMT Python package for creating Surrogate models

https://smt.readthedocs.io/en/latest/_src_docs/surrogate_models.html


Try our SMT Python package for doing Efficient Global Optimization 

https://smt.readthedocs.io/en/latest/_src_docs/applications/ego.html

****************

SVD

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




