function xi = poles(obj)
%POLES    Return the poles of an RKFUN.
  
  xi = util_pencil_poles(obj.K, obj.H);
  xi = xi(:);
  
  % Only output the m smallest poles, where m is denominator
  % degree.
  t = type(obj);
  [~, ind] = sort(abs(xi), 'ascend');
  xi = xi(ind);
  xi = xi(1:t(2));
  
end
