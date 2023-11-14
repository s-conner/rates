/*
Example analysis
Marginal versus conditional rate estimation for count and recurrent event data with an estimand framework
Sarah Conner, Yijie Zhou, and Tu Xu
November 13 2023
*/

proc import out=dat datafile='E:\users\conner\PEx rate\example.csv' dbms=csv replace; getnames=yes; run;

data dat; 
set dat;
log_t=log(t);
run;


* Poisson;
proc genmod data=dat;
class trt(ref='0');
model y=trt x / dist=poisson link=log offset=log_t;
lsmeans trt / cl diff exp;
run;

* NegBin;
proc genmod data=dat;
class trt(ref='0');
model y=trt x / dist=nb link=log offset=log_t;
lsmeans trt / cl diff exp;
run;
