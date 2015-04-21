# CCBR Git Tutorial

This document originated with a presentation given at the CCBR Tech-Dev meeting on February 25, 2015. Most of the information in this document (and a lot more) can be found at https://www.atlassian.com/git/tutorials/.

## Retrieving information from a remote repository
`git clone <location of repository to clone> <location to write the cloned repository>` clones a git repository to the location indicated. The location must either be an empty directory or must not exist -- git will create the needed directories

`git fetch` allows you to review what has been done on a remote branch without having to integrate the changes into your local repository.
These changes can be reviewed using git log (see below) and then merged into your local directory using git merge. git pull essentially rolls these two steps into one, fetching and merging all updates into your local copy of the branch.

`git pull` updates your local copy of the repository with the latest remote version.
This will fail if there are files with changes in both the local and remote repositories.

`git pull --rebase` pull all new commits and merges them into your local version, but rebases any changes on top of what anyone else has done (i.e. this isn’t really a merge).

`git checkout <commit> <file (optional)>` - View an old state of the repository without altering the local state. This will temporarily alter the files you will see in the local repository.
Checking out a particular file will modify the current state of your local repository, and it will appear in the list of files that have been changed and are awaiting committal.

`git checkout <branch>` will let you switch to a different branch.
Any development on a CCBR pipeline repository should be done in a branch other than the main branch, preferably on a fork of the main repository. The -b option can be added to create the branch indicated.

`git checkout master` will restore your local repository to the current production version.

`git remote -v` shows all remote repositories that are connected with your local repository

`git remote add <name> <url>` allows you to add a link to another remote repository located at <url> and named <name>.

Restoring information when you want to reset a local repository or revert an earlier commit. This sequence can be used to reset (and delete!) all changes you have made. If you really mess things up on your local copy, you can use this to get a fresh start.
```
git fetch
git reset --hard origin/master
```

`git revert <commit>` undoes a previous commit. This actually creates a new commit with the changes from the relevant commit redacted, thus preserving the history of the repository. 
For a figure demonstrating the difference between resetting and reverting, [click here](https://www.google.com/url?q=https%3A%2F%2Fwww.atlassian.com%2Fgit%2Ftutorials%2Fmerging-vs-rebasing%2Fconceptual-overview&sa=D&sntz=1&usg=AFQjCNHrmVA-uJGYP0W53fDwcPrlHJnW7Q).

`git clean -n` does a dry run of git clean, showing what files will be removed by git clean.

`git clean -f` removes all untracked files in your local repository.

It is possible to rewrite the history of your repository, but this shouldn’t be something we often come across and is beyond the scope of this document (except perhaps when rebasing a feature branch as discussed below). For more information on this [click here](https://www.google.com/url?q=https%3A%2F%2Fwww.atlassian.com%2Fgit%2Ftutorials%2Frewriting-history&sa=D&sntz=1&usg=AFQjCNG9vA_tUSCzabbGGVU1SeGFzPMlfQ).

## Inspecting a repository
`git branch` lists all branches in the repository

`git status` gives you a summary of the untracked files, changes you have made, those that have been added, and those that have been committed (but not pushed).
This will also show you conflicts when trying to merge two branches.

`git log` displays the entire commit history.

`git log -n 10` displays first 10 lines of the commit history

`git log --oneline` displays a condensed one-line message for each commit.

`git log -p` displays the full diff information for each commit (verbose output). For more information on git log see [this article](https://www.atlassian.com/git/tutorials/git-log/).

## Adding your modifications to a remote repository
Once you have made changes to files that you want to commit back to the source repository, add, commit and push the files. You can add some or all changed files to a commit.
Upon using git commit, you will be prompted to enter a short description of your commit. The general best practice is to enter a short (<50 characters) message, followed by a blank line and more detail if desired.
Multiple commits may be made before pushing, but changes will not be available to other developers until they are pushed to the source repository.
All files and directories represented in .gitignore will not be tracked or pushed to the repository.
```
git add <files or directories to add>
git add --all
git commit
git push
```

## Creating and managing branches
`git branch <name of new branch>` creates a new branch from the current commit. This can be done from a previous commit if desired by checking out the previous commit and branching from there.
The new branch is not checked out.

`git branch -b <name of new branch>` will create a new branch and check it out (this doesn’t seem to be a valid option in OS X’s implementation need to double check this syntax, but `git checkout -b` works).

`git merge <branch>` will merge <branch> with the current branch.
A fast-forward merge takes place when the root of one branch is at the head of the other branch (see [this image](https://www.google.com/url?q=https%3A%2F%2Fwww.atlassian.com%2Fgit%2Fimages%2Ftutorials%2Fcollaborating%2Fusing-branches%2F07.svg&sa=D&sntz=1&usg=AFQjCNHiGdgEC3pRL58UVa6u6ELSy2Z0eA)).
A 3-way merge takes place when changes have been made to both branches. Conflicts may arise. Use git status to identify where the conflicts are and manually fix them (see [this image](https://www.google.com/url?q=https%3A%2F%2Fwww.atlassian.com%2Fgit%2Fimages%2Ftutorials%2Fcollaborating%2Fusing-branches%2F08.svg&sa=D&sntz=1&usg=AFQjCNENQFRniNkM6Dsz3Jx7vKvcuP8x8g)).

Using `git rebase` is a great way to keep your branch up to date with the master branch. It effectively re-roots your branch at the head of the master branch. This can also be accomplished with git merge, but it complicates the history of your branch (there are pos and cons to each approach -- see [this article](https://www.google.com/url?q=https%3A%2F%2Fwww.atlassian.com%2Fgit%2Ftutorials%2Fmerging-vs-rebasing%2Fconceptual-overview&sa=D&sntz=1&usg=AFQjCNHrmVA-uJGYP0W53fDwcPrlHJnW7Q) for more details):
```
git checkout <branch>
git rebase master
```

`git stash` stashes changes when you aren’t ready to commit when you need to do some work on a different branch -- This is very useful if you have a couple of different active branches you are working on, or if you would like to checkout a different branch before you are ready to commit your changes to the current branch you are working on (see [this article](http://www.google.com/url?q=http%3A%2F%2Fgit-scm.com%2Fbook%2Fen%2Fv1%2FGit-Tools-Stashing&sa=D&sntz=1&usg=AFQjCNH_vKrIFwUO9iJM1yuinfsk2HW5eQ) for more details and examples). This stashes all uncommitted work, including changes you have staged using git add.

`git stash list` will show you all active stashes you have on the stack.

`git stash apply <stash - optional>` will apply a stash back to your workspace to allow you to continue where you left off. By default it will automatically apply the most recent stash in the stack, but <stash> can be supplied if you want to apply a different stash.

## SourceTree
[Sourcetree](http://sourcetreeapp.com) is a nice piece of software that provides a GUI interface for administering git repositories.
