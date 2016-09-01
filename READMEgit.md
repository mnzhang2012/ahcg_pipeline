# Readme file for change the remote url for git repositories on Mac terminal
1. Change the current working directory to local git repository folder
2. List exiting remotes for the repository
   $ git remote -v
3. Use the git remote set-url command to change the remote repository's URL
   $ git remote set-url origin https://github.com/mnzhang2012/ahcg_pipeline.git
4. Verify the remote URL has changed
   $ git remote -v