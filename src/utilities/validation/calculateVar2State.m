function [r, v_mRNA, v_ML] = calculateVar2State(kon,koff,c_true)

  r = c_true*kon / (c_true*kon + koff);      
  v_mRNA = 2*(c_true*kon+koff) / (c_true*kon*koff);
  v_ML = (c_true*kon+koff)/(c_true*kon*koff);