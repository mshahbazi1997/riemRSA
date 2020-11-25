function COV = covPattern(Pattern, n_RIO, ROIs_load_addr)

returnHear = pwd;
cd(['..' filesep '..'])
ROIs = niftiread(ROIs_load_addr);
cd(returnHear)

Pattern = Pattern(ROIs(:) == n_RIO, :);
COV = cov(Pattern);

end
