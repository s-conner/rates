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


* Poisson w/ G-estimation;
proc genmod data=dat;
class trt(ref='0');
model y=trt x / dist=poisson link=log offset=log_t;
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
var pred1;
output out=rate1;
run;
proc means data=preds0 mean;
var pred0;
output out=rate0;
run;

proc sql;
select pred1 as rate1, pred0 as rate0, pred1/pred0 as rr
from rate1 r1
join rate0 r0 on 1
where r1._STAT_='MEAN' and r0._STAT_='MEAN';
quit;

* Bootstrap 95% CI;
