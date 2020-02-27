 function [av]=isoutlier(FF)
  mk=median(FF);
  M_d=mad(FF,1);
  c=-1/(sqrt(2)*erfcinv(3/2));
  smad=c*M_d;
  tsmad=3*smad;
  av=(abs(FF-mk)>=tsmad);
 end