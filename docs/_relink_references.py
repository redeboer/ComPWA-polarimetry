# cspell:ignore refdomain refspecific reftarget reftype
"""Abbreviated the annotations generated by sphinx-autodoc.

It's not necessary to generate the full path of type hints, because they are
rendered as clickable links.

See also https://github.com/sphinx-doc/sphinx/issues/5868.
"""
from __future__ import annotations

import sphinx.domains.python
from docutils import nodes
from sphinx.addnodes import pending_xref, pending_xref_condition
from sphinx.domains.python import parse_reftarget
from sphinx.environment import BuildEnvironment

__TARGET_SUBSTITUTIONS = {
    "DataSample": "tensorwaves.interface.DataSample",
    "Literal[(-1, 1)]": "typing.Literal",
    "Literal[- 1, 1]": "typing.Literal",
    "OuterStates": "polarization.decay.OuterStates",
    "ParametrizedBackendFunction": "tensorwaves.function.ParametrizedBackendFunction",
    "Path": "pathlib.Path",
    "Pattern": "typing.Pattern",
    "PoolSum": "ampform.sympy.PoolSum",
    "UnevaluatedExpression": "ampform.sympy.UnevaluatedExpression",
    "implement_doit_method": "ampform.sympy.implement_doit_method",
    "sp.Expr": "sympy.core.expr.Expr",
    "sp.Indexed": "sympy.tensor.indexed.Indexed",
    "sp.Rational": "sympy.core.numbers.Rational",
    "sp.Symbol": "sympy.core.symbol.Symbol",
    "sp.acos": "sympy.functions.elementary.trigonometric.acos",
}
__REF_TYPE_SUBSTITUTIONS = {
    "polarization.decay.OuterStates": "obj",
    "tensorwaves.interface.DataSample": "obj",
}


def _new_type_to_xref(
    target: str,
    env: BuildEnvironment = None,
    suppress_prefix: bool = False,
) -> pending_xref:
    reftype, target, title, refspecific = parse_reftarget(target, suppress_prefix)
    target = __TARGET_SUBSTITUTIONS.get(target, target)
    reftype = __REF_TYPE_SUBSTITUTIONS.get(target, reftype)
    assert env is not None
    return pending_xref(
        "",
        *__create_nodes(env, title),
        refdomain="py",
        reftype=reftype,
        reftarget=target,
        refspecific=refspecific,
        **__get_env_kwargs(env),
    )


def __get_env_kwargs(env: BuildEnvironment) -> dict:
    if env:
        return {
            "py:module": env.ref_context.get("py:module"),
            "py:class": env.ref_context.get("py:class"),
        }
    return {}


def __create_nodes(env: BuildEnvironment, title: str) -> list[nodes.Node]:
    short_name = title.split(".")[-1]
    if env.config.python_use_unqualified_type_names:
        return [
            pending_xref_condition("", short_name, condition="resolved"),
            pending_xref_condition("", title, condition="*"),
        ]
    return [nodes.Text(short_name)]


def relink_references() -> None:
    sphinx.domains.python.type_to_xref = _new_type_to_xref
