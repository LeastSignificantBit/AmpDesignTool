function s=num2eng(d, prec=4)
  if d>=1e9
    s=[num2str(d/1e9,prec) ' G'];
  elseif d>=1e6
    s=[num2str(d/1e6,prec) ' M'];
  elseif d>=1e3
    s=[num2str(d/1e3,prec) ' k'];
  elseif d>=1
    s=[num2str(d,prec) ' '];
  elseif d>=1e-3
    s=[num2str(d/1e-3,prec) ' m'];
  elseif d>=1e-6
    s=[num2str(d/1e-6,prec) ' µ'];
  elseif d>=1e-9
    s=[num2str(d/1e-9,prec) ' n'];
  elseif d>=1e-12
    s=[num2str(d/1e-12,prec) ' p'];
  elseif d>=1e-15
    s=[num2str(d/1e-15,prec) ' f'];
  else
    s=[num2str(d/1e-18,prec) ' a'];
  endif
endfunction