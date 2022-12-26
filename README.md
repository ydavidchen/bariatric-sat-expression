# Gene expression analysis of subcutaneous adipose tissue (SAT) exposed to bariatric surgery

Status: Peer-reviewed and revised

Author: David Chen, Ph.D.

## Introduction

Bariatric surgery is currently the most effective approach to achieve weight loss as an obesity intervention. However, the molecular mechanisms behind bariatric surgery remains elusive, especially given that the existing studies tend to have very small sample sizes. 

## Design

Five existing datasets from Gene Expression Omnibus were re-processed, joined, and variance-filtered. The resulting dataset has 5,000 genes across 239 SATs (126 post-surgery, 111 pre-surgery/controls).

Unsupervised and supervised analyses surrounding pre/post-bariatric surgery status were then performed.

## Source Code

End-to-end R code is available here.

* Individual dataset cleaning
* Combined data processing, normalization, and batch-effect removal
* Unsupervised analyses
* Supervised analysis

Downstream biological interpretations (Gene Ontology and KEGG) are performed on the [WebGestalt server](www.webgestalt.org).

## Supplemental Results

The complete list of genes tested for differential expression can be found in _supplemental results_ subdirectory.

Copy right &copy; Y. David Chen. All rights reserved.

Citation:



