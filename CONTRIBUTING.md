Contributing to graph-imputation-match
===

1) Log into Github web interface with your username-nmdp account

2) Browse to the repo at https://github.com/nmdp-bioinformatics/graph-imputation-match, hit the Fork button

3) Copy the clone URL from the Github web page for the fork (something like https://github.com/username-nmdp/graph-imputation-match.git)

4) Clone the fork
```bash
$ git clone https://github.com/username-nmdp/graph-imputation-match.git
$ cd graph-imputation-match
```

5) Add upstream as remote
```bash
$ git remote add upstream https://github.com/nmdp-bioinformatics/graph-imputation-match
```

6) Pull and merge latest changes from upstream master to your local master branch
```bash
$ git checkout master
$ git pull upstream master
$ git push
```

7) Create a new local feature branch
```bash
$ git checkout -b new-feature-branch
```

8) Edit files locally

9) Commit changes to local feature branch
```bash
$ git commit -m "made changes"
```

10) Push changes from local feature branch to remote feature branch on your fork
```bash
$ git push origin new-feature-branch
```

11) Browse to the Github web page for your fork repo (something like https://github.com/username-nmdp/graph-imputation-match) and hit the new pull request button

12) Edit the pull request description and hit save (or create or new, I can't remember) button

13) Other contributors will review the changes in the pull request

14) When the pull request looks good, it is merged into the upstream repo

15) Hit the delete branch button to delete your remote feature branch (the commits have been merge upstream, so it is no longer necessary)

16) Delete your local feature branch
```bash
$ git branch -d new-feature-branch
```
