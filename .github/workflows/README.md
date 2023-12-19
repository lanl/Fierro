# Github Actions

Github actions can do a lot of different things. There is even an opensource Github action [marketplace](https://github.com/marketplace?category=&query=&type=actions&verification=) that hosts tons of actions that you can easily add to a project. However, we currently only use Github actions for one thing - publishing Anaconda packages.

## Action Overview
Github actions are processes that are controlled by yaml configuration scripts under `/.github/workflows` (here). Each file in this folder defines an action, and each action typically runs independently from one another[^1]. An action is made of two primary components: the triggers and the jobs.


[^1]: There are a couple of ways to chain actions to split things like building and uploading into separate processes with a dependency.

### Triggers
Actions automatically run when one of their [triggering events](https://docs.github.com/en/actions/using-workflows/events-that-trigger-workflows) occurs. For us, the important ones here on push, and on workflow_dispatch. By specifying "on.push" and "on.workflow_dispatch" in the yaml file, we can tell github to run this action when either a push occurs or we manually dispatch the workflow from a button in the github actions webview.

```yaml
on: 
  push:
    paths:
      - .conda/voxelizer/**
      - src/Voxelizer/**
      - python/Voxelizer/**
      - .github/workflows/publish-voxelizer.yaml
  workflow_dispatch:
```

The above "on" block tells github to trigger the action when a push occurs that matches the listed file patterns. In this case, if the Voxelizer source code changes, or the anaconda build metadata changes for the package.

### Jobs
A job in a github action can do almost anything, and consists of two primary components: the run strategy and the run steps. Just like in conda-build, the run strategy can define a matrix of options for which to run the steps. In our action designed to build and publish anaconda packages, we have specified the following matrix to tell github to run this action on both MacOS and Linux.
```yaml
jobs:
    ...
    strategy:
      matrix:
        include:
          - os: macos-latest
          - os: ubuntu-latest
    runs-on: ${{ matrix.os }}
    ...
```

### Secrets
Lastly, it is often necessary to give a github runner (process that executes actions) special permissions or access to restricted resources. In our case, it is an Anaconda channel. We want the runner to publish new builds of Anaconda packages. This is typically done by giving the runner a token that has restricted access to an account to do only what it needs to do. However, we don't want to release the token to the public, so instead we use a repository secret. These are hidden variables that can only be read by github runners (not even humans can see them), but can be set by an authorized user. We use `ANACONDA_TOKEN` as a limited access token that the runner can use by putting the correct template into the action configuration file: `${{ secrets.ANACONDA_TOKEN }}`. That will get replaced by the secret during execution.

## Our Actions
Since we only do one kind of action (publishing Anaconda packages), we only need 2 kinds of action files. The first (`/.github/workflows/build-conda-package.yaml`) is a reusable action that doesn't automatically trigger on its own, but instead is invoked by other actions. This action takes two configurable arguments: 
```yaml
on: 
  workflow_call:
    inputs:
      recipe_dir:
        required: true
        type: string
      additional_variant:
        required: false
        type: string
        default: '.conda/empty.yaml'
```

The recipe_dir is the folder containing the conda build `meta.yaml` file as well as license information. The additional_variants is a file path to another `conda_build_config.yaml` file to add into the build matrix. As mentioned in `.conda/README.md`, we use this to add specific python version variants for the voxelizer package.

The second is the package specific action that triggers off of package specific source code changes and invokes `build-conda-package` with the correct package recipe location and additional build variants. In the case of our source-in-repo packages, these trigger of off relevant source code changes. For other dependencies, where the source is not contained in the repo, we configure the actions for manual dispatch.


# Debugging
Complex github actions can be very annoying to debug for two reasons:
1. The feedback loop is long/slow. Pushing changes to github and waiting for a runner to initialize can easily add ~5 mins to your feedback loop.
2. You don't have easy access to the runtime environment. Github automatically will unmount the docker image and be on its way if an action fails.

We have some tools to deal with this, though. If you can run docker (not WSL1), you can use [act](https://github.com/nektos/act) to execute github actions locally. This is **much** faster for debugging things like action syntax. We also have a tool to deal with #2 that I have built into the `build-conda-package` action. We have a `debug` option that enables a special workflow step after the package building step. If the debug flag is set, we make use of the [tmate](https://github.com/mxschmitt/action-tmate) github action to create a temporary ssh host and wait for 30 minutes. On the github Actions page you can inspect the log and find something like:
```
SSH: ssh bJV7VnwDn5NUzBn6fGrkJYQg9@nyc1.tmate.io
```

If you are the user that triggered the action, you can ssh into the docker container that was runing your build action. From there you can continue debugging.