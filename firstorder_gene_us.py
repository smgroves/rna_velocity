# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 18:27:53 2019

@author: johnv
"""

from pysb import Model, Monomer, Parameter, Rule, Initial, Observable

Model()


# ------------------------------------------------------
# DEFINE MONOMERS

Monomer('g', ['b'])
Monomer('u')
Monomer('s')
Monomer('p', ['b'])

# ---------------------------------------------------------
# Define parameters

Parameter('k_0', 1.0)     # transcription rates
Parameter('k_1', 40.0)

Parameter('k_p', 3.0)     # translation rate

Parameter('beta', 2.0)    # unspliced to spliced rate

Parameter('d_s', 1.0)     # degradation rates
Parameter('d_p', 1.0) 

Parameter('b_f', 0.1)     # gene-protein forward/backward binding rates
Parameter('b_r', 10.0)

# ------------------------------------------------

Rule('transcription_0bound', g(b=None) >> g(b=None) + u(), k_0)              # Transcription reactions
Rule('transcription_1bound', g(b=1) % p(b=1) >> g(b=1) % p(b=1) + u(), k_1)

Rule('translation', s() >> s() + p(b=None), k_p)    # Translation reaction

Rule('splicing', u() >> s(), beta)        # splicing reaction

Rule('death_s', s() >> None, d_s)                 # Degradation reactions
Rule('death_p', p(b=None) >> None, d_p)

Rule('g_bind', g(b=None) + p(b=None) | g(b=1) % p(b=1), b_f, b_r)      # gene-protein binding reaction


# ------------------------------------------------


# initial conditions
Parameter('g_0',    1)
Parameter('u_0',  0)
Parameter('s_0', 0)
Parameter('p_0',  0)

Initial(g(b=None), g_0)
Initial(u, u_0)
Initial(s, s_0)
Initial(p(b=None), p_0)


# Observables
Observable('o_g_0', g(b=None))
Observable('o_g_1', g(b=1) % p(b=1))
Observable('o_u', u)
Observable('o_s', s)
Observable('o_p', p(b=None))