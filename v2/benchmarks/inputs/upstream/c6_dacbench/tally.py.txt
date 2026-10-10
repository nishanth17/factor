#!/usr/bin/env python3

nodes = 0
nodesdetail = {}

def clear():
  global nodes
  global nodesdetail
  nodes = 0
  nodesdetail = {}

def node(what):
  global nodes
  nodes += 1
  if what not in nodesdetail:
    nodesdetail[what] = 0
  nodesdetail[what] += 1
