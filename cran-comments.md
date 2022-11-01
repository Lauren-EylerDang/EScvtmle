For a CRAN submission we recommend that you fix all NOTEs, WARNINGs and ERRORs.
## Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)

## R CMD check results
❯ On windows-x86_64-devel (r-devel)
  checking CRAN incoming feasibility ... [14s] NOTE
  Maintainer: 'Lauren Eyler Dang <lauren.eyler@berkeley.edu>'
  
  New submission
  
  Possibly misspelled words in DESCRIPTION:
    CVTMLE (11:96, 11:385)
    RCT (4:9, 11:302, 11:319, 11:622)
    TMLE (3:31)
    al (11:905)
    et (11:902)
    
Note: This is a new submission. The possibly misspelled words are either common (e.g. et al. ) or previously defined abbreviations.

❯ On windows-x86_64-devel (r-devel)
  checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

❯ On ubuntu-gcc-release (r-release)
  checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Lauren Eyler Dang <lauren.eyler@berkeley.edu>’
  
  New submission
  
  Possibly misspelled words in DESCRIPTION:
    al (11:905)
    CVTMLE (11:96, 11:385)
    et (11:902)
    RCT (4:9, 11:302, 11:319, 11:622)
    TMLE (3:31)
    
Note: This is a new submission. The possibly misspelled words are either common (e.g. et al. ) or previously defined abbreviations.

❯ On fedora-clang-devel (r-devel)
  checking CRAN incoming feasibility ... [6s/25s] NOTE
  Maintainer: ‘Lauren Eyler Dang <lauren.eyler@berkeley.edu>’
  
  New submission
  
  Possibly misspelled words in DESCRIPTION:
    CVTMLE (11:96, 11:385)
    RCT (4:9, 11:302, 11:319, 11:622)
    TMLE (3:31)
    al (11:905)
    et (11:902)
    
Note: This is a new submission. The possibly misspelled words are either common (e.g. et al. ) or previously defined abbreviations.

❯ On fedora-clang-devel (r-devel)
  checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found

0 errors ✔ | 0 warnings ✔ | 5 notes ✖
