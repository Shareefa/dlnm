# Research project on the benefits of an optimized Distributed Lag Non-linear Model

## Problem Statement

Does an optimized DLMN outperform standard DLMN with mixmeta on large datasets?

## Background

mixmeta repo: https://github.com/gasparrini/mixmeta
mixmeta paper: http://www.ag-myresearch.com/uploads/1/3/8/6/13864925/2019_sera_statmed.pdf

This is a repo for a DLMN model for epidemiology studies. This implementation is limited a few ways details in `benchmarks/MODEL_EXPLAINER.md`.

Mixmeta is a solution for the dlmn dataset limitations. 



## Goals

* Plan list of optimizations for DLNM
* Implement optimizations
* Benchmark the performance gains of optimizing DLNM vs unoptimized DLNM
* Generate synthetic >10GB data to use with mixmeta and optimized DLNM
* Benchmark the performance gains of optimized DLMN vs unoptimized DLNM with mixmeta along a few different configurations

## Output

* Generate a report on how well the model converged on the underlying data under different configs.
