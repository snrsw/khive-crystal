# khive-crystal

This repository contains the experimetns code for our paper `Algorithms for crystal structure on K-hives of type A`.

## Overview

Let $\mathfrak{g}$ be a Lie algebra of type $A_{n-1}$. Let $P^+$ be a dominant integral weight lattice of $\mathfrak{g}$. For $\lambda \in P^+$, let $\mathbb{H}(\lambda)$ be a set of K-hives which right edge labels is determied by $\lambda$.
`khive-crystal` compute the crystal structure on $\mathbb{H}(\lambda)$.

## Getting started

### Install with poetry

```bash
pip install git+https://github.com/snrsw/khive-crystal@master
```

## Usage

For more information, see https://snrsw.github.io/khive-crystal/

### K-hive

```python
>>> from khive_crystal import khive
>>> khive(alpha=[2, 1, 0], beta=[2, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
KHive(alpha=[2, 1, 0], beta=[2, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])
```

```python
>>> from khive_crystal import khive, view
>>> view(khive(alpha=[2, 1, 0], beta=[2, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]]))
```

![](<.github/images/KHive(n=3,%20alpha=[2,%201,%200],%20beta=[2,%201,%200],%20gamma=[0,%200,%200],%20Uij=[[0,%200],%20[0]]).png>)

### Crystal structure

#### $\mathbb{H}(\Lambda_k)$

For $k = 1, 2, \dots, n-1$, let $\Lambda_k$ be the fundamental weight of $\mathfrak{g}$.
The crystal structure on $\mathbb{H}(\Lambda_k)$ is computed as follows.

```python
>>> from khive_crystal.khive import KHive
>>> from khive_crystal import khive, view, f, phi, e, epsilon
>>> H: KHive = khive(
...     n=3,
...     alpha=[1, 1, 0],
...     beta=[1, 1, 0],
...     gamma=[0, 0, 0],
...     Uij=[[0, 0], [0]]
... ) # alpha = \Lambda_2
>>> phi(i=2)(H)
1
>>> f(i=2)(H)
KHive(alpha=[1, 1, 0], beta=[1, 0, 1], gamma=[0, 0, 0], Uij=[[0, 0], [1]])
>>> epsilon(i=2)(H)
0
>>> e(i=2)(H)
# None
```

#### $\mathbb{H}(\lambda)$

Let $\lambda \in P^+$.
The crystal structure on $\mathbb{H}(\lambda)$ is computed in two ways.

##### As a submodule of tensor products of $\mathbb{H}(\Lambda_k)$

The crystal structure on $\mathbb{H}(\lambda)$ is determined such that the map $\Psi \colon \mathbb{H}(\lambda) \to \bigotimes_{k} \mathbb{H}(\Lambda_k)$ is a crystal embedding.

```python
>>> from khive_crystal.khive import KHive
>>> from khive_crystal import khive, view, f, phi, e, epsilon, psi, psi_lambda, psi_inv
>>> H: KHive = khive(n=3, alpha=[3, 1, 0], beta=[2, 2, 0], gamma=[0, 0, 0], Uij=[[1, 0], [0]])
>>> psi_lambda(H=H)
[KHive(alpha=[2, 0, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[1, 0], [0]]), KHive(alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])]
>>> psi(H=H)
[KHive(alpha=[1, 0, 0], beta=[0, 1, 0], gamma=[0, 0, 0], Uij=[[1, 0], [0]]), KHive(alpha=[1, 0, 0], beta=[1, 0, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]]), KHive(alpha=[1, 1, 0], beta=[1, 1, 0], gamma=[0, 0, 0], Uij=[[0, 0], [0]])]
>>> psi_inv(
...     H=f(i=1)(
...         psi(
...             H=H
...         )
...     )
... )
KHive(alpha=[3, 1, 0], beta=[1, 3, 0], gamma=[0, 0, 0], Uij=[[2, 0], [0]])
```

##### Explicit method

```python
>>> from khive_crystal import khive, view, f, phi, e, epsilon
>>> H: KHive = khive(n=3, alpha=[3, 1, 0], beta=[2, 2, 0], gamma=[0, 0, 0], Uij=[[1, 0], [0]])
>>> f(i=2)(H)
KHive(alpha=[3, 1, 0], beta=[1, 3, 0], gamma=[0, 0, 0], Uij=[[2, 0], [0]])
```

#### Crystal graph

```python
>>> from khive_crystal import khives, crystal_graph
>>> crystal_graph(khives(n=3, alpha=[2, 1, 0]))
```

![](<.github/images/khives(n=3,%20alpha=[2,%201,%200]).png>)
