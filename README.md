# Reproducible paper proof of concept

## Purpose

To test approaches for building a completely reproducible document. This repository should act as something of a sandbox to experiment with methods and to practice using tools.

## Strategy

To be successful, the document should be able to be compiled in a fresh environment by cloning this repository and running a script/makefile.

This will require reproducible:
- data access
- compute environment
- code/scripts
- document generation

### Tools

As a starting point, the following tools will be used:

- Git. To manage the contributions of multiple people we will attempt to follow this model (http://nvie.com/posts/a-successful-git-branching-model/) for branching. tldr: don't make commits on the master branch. master only gets modified via pull requests.

- [Bookdown](https://bookdown.org/yihui/bookdown/) - an R package that enables writing technical documents in markdown

- [Singularity](https://github.com/singularityware/singularity) - a container system for reproducible/portable compute environments
>>>>>>> readme

