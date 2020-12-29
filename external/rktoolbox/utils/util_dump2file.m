function util_dump2file(filename, varname, var)
% UTIL_DUMP2FILE    Appends data into a *.mat file.
%
% util_dump2file(filename, varname, var) stores the data var into
% the file filename using the name varname within the file. This
% file is mainly used to export data to Python for plotting
% purposes.  

  tmp = struct(varname, var);
  try
    save(filename, '-struct', 'tmp', '-append');
  catch
    save(filename, '-struct', 'tmp');
  end
  
end