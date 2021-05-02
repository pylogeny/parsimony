# Parsimony Phylogenetics with Python

This package offers a simple parsimony implementation with step matrices.

## Installation

```
$ pip install pyloparsimony
```

## Usage

```python
>>> from pyloparsimony import parsimony, EXAMPLES, print_scenario
>>> tree_weights, scenarios, weights = parsimony(EXAMPLES["e1"])
>>> print_scenario(scenarios[0], tree=EXAMPLES["e1"]["tree"])
                              ┌─A/b
                    ┌─Edge1/b─┤
                    │         └─B/c
          ┌─Edge3/b─┤
          │         │         ┌─C/a
──Root/b──┤         └─Edge2/a─┤
          │                   └─D/a
          └─E/b

```
