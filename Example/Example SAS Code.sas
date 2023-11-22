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

t2=rand('UNIFORM',0.8,1.2);
log_t2=log(t2);
run;


* Poisson;
proc genmod data=dat;
class trt(ref='0');
model y=trt x / dist=poisson link=log offset=log_t2;
lsmeans trt / cl diff exp;
run;

* NegBin;
proc genmod data=dat;
class trt(ref='0');
model y=trt x / dist=nb link=log offset=log_t2;
lsmeans trt / cl diff exp;
run;


* Poisson w/ G-estimation;
proc genmod data=dat;
class trt(ref='0');
model y=trt x / dist=poisson link=log offset=log_t2;
lsmeans trt / cl diff exp;
store out=mod;
run;

data trt1; set dat; trt=1; run;
data trt0; set dat; trt=0; run;

proc plm source=mod;
score data=trt1 out=preds1 pred=pred1 / ilink;
run;
proc plm source=mod;
score data=trt0 out=preds0 pred=pred0 / ilink;
run;

proc means data=preds1 mean;
var pred1 t2;
output out=rate1 mean=pred1 t2;
run;
proc means data=preds0 mean;
var pred0 t2;
output out=rate0 mean=pred0 t2;
run;

proc sql;
select r1.pred1/r1.t2 as rate1, r0.pred0/r0.t2 as rate0, r1.pred1/r0.pred0 as rr
from rate1 r1
join rate0 r0 on 1;
quit;

* Bootstrap 95% CI;
