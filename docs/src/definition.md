# Definition

We first introduce a base-`B` representation with `B = 2, 3, 4, ...`.
Throughout this guide, quantics digits and grid indices are 1-based.

We represent a positive integer `X >= 1` as

$$
X = \sum_{i=1}^{R} (x_i - 1) B^{R-i} + 1,
$$

where each digit `x_i` satisfies `1 <= x_i <= B` and `R` is the number of
digits. In this crate, the base-`B` representation of `X` is stored as the
vector

$$
[x_1, \ldots, x_R].
$$

For multiple variables, the crate supports fused and interleaved unfolding
schemes. For three variables `X`, `Y`, and `Z`, suppose their base-`B`
representations are

$$
[x_1, \ldots, x_R], \quad [y_1, \ldots, y_R], \quad [z_1, \ldots, z_R].
$$

The interleaved representation is

$$
[x_1, y_1, z_1, x_2, y_2, z_2, \ldots, x_R, y_R, z_R].
$$

The fused representation is

$$
[\alpha_1, \alpha_2, \ldots, \alpha_R],
$$

where

$$
\alpha_i = (x_i - 1) + B (y_i - 1) + B^2 (z_i - 1) + 1
$$

and therefore

$$
1 \le \alpha_i \le B^3.
$$

In fused ordering, the `x` digit runs fastest at each digit level. This matches
the convention used by the Julia package and generalizes to any number of
variables.
