%%
% Runs all the unit tests in the 'tests' directory.

rng('default')

fprintf('\n Rational Krylov Toolbox version 2.4.\n Please wait while running all tests.\n')
fprintf(' %s\n', repmat('-', 1, 36))
d = dir('*');
folder = [d(:).isdir];
folder = {d(folder).name}';
folder(ismember(folder, {'..'})) = [];

tcount = 0;
failed = 0;
tic
for fid = 1:length(folder)
  cd(folder{fid})
  tests = dir('test_*.m');
  if length(tests) == 0, continue; end
  fprintf(' Folder: %s\n', folder{fid})  
  for j = 1:length(tests)
    tcount = tcount +1;
    tid = sprintf('[%d]', tcount);

    fn = tests(j).name; 
    fn = fn(1:end-2);        
    
    fprintf('  %5s %s %s ', ...
            tid, ...
            fn(6:end), ...
            repmat('.', 25-length(fn), 1))
    
    try,     check = feval(fn);
    catch e, check = 0;         end
    
    if ~all(check), fprintf('failed. <---\n'), failed = failed + 1;      
    else            fprintf('passed.\n'),      end
  end
  if fid ~= length(folder), fprintf('\n'), end
  if ~strcmp(folder{fid}, '.'), cd('..'), end
end % fid = 1:length(folder)
t = toc;
fprintf(' %s\n', repmat('-', 1, 36))

if failed, fprintf(' %2d tests failed. \n\n', failed)
else       fprintf(' All tests passed in %2.2f seconds.\n\n', t);
end
