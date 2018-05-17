---
title       : Propensity score calipers and the overlap condition
author      : Ben Hansen, UMich Statistics (bbh@umich.edu).
date        : IISA, May 2018
job         : benthestatistician.github.io/PISE (slides)
framework   : io2012        # {io2012, html5slides, shower, dzslides, ...}
highlighter : highlight.js  # {highlight.js, prettify, highlight}
hitheme     : tomorrow      #
widgets     : [mathjax]            # {mathjax, quiz, bootstrap}
mode        : selfcontained # {standalone, draft}
--- 

# Outline

1. <span class = 'red'>Problem: determining the region of  overlap</span>
2. Propensity score caliper matching
3. PPSE calipers
4. Applications & discussion    




--- &twocol    

## Overlap & bias in comparisons w/ 1 covariate

*** =left

> - The quasi-experimental setup: $X$, $Y$; $Z \in \{0,1\}$, $Y_c$.
> - <embed height="400px" width="400px" src="./images/pretest-comp-group-design.jpg">

*** =right

> - Overlap? 
>     - Observations well outside region of overlap should be excluded.
>     - Wasteful to exclude members of $\{ i: Z_i=1\}$ because  their $x$ is just outside $\mathrm{range}(\{x_j: Z_j=0\})$.
>     -  May distort research question, too.
>     - Close/far distinctions presuppose a model.

> - Balance - i.e., is $\bar{x}_t \approx \bar{x}_c$?
>    - The question presupposes overlap.  (OW mean comparisons can be misleading.)
>    - A topic for another day...

---

## Statistical assumptions for general $X$es

>- Now, multiple $X$es: $\mathbf{X}$, not $X$.
>- Strong ignorability conditions:
>    - No unmeasured confounding: $Y_c \perp Z | \mathbf{X}$.
>    - Overlap: $\mathrm{Pr}(Z=1  | \mathbf{X}) < 1$.      
>- No unmeasured confounding says, there are natural experiments in each "cell" of $\mathbf{X}$.
>- Problem: if multiple ${X}$es, many observations will be all alone in their cells.
>- Rosenbaum & Rubin's (1983) solution is match on $\widehat{\mathrm{Pr}}(Z=1| \mathbf{X}=\mathbf{x})$, the estimated propensity score (PS).
>- Second problem: In a cell $\mathbf{x}_0$ s.t. $\mathrm{Pr}(Z=1| \mathbf{X}=\mathbf{x}_0)=1$, there's no data to estimate $\mathrm{Pr}(Y_c \in \cdot |  \mathbf{X}=\mathbf{x}_0)$, even as $n \uparrow \infty$.  (Whether you match on the PS, weight for it, ...)
>- The second problem motivates the method of strict overlap, i.e. reducing analytic sample to region of common support on $\widehat{\mathrm{PS}}$: $$\cap_{z=0,1}[\min(\{\widehat{PS}_i:Z_i=z\}), \max(\{\widehat{PS}_i:Z_i=z\})] .$$ 


--- &twocol
##  Example 1: violence & public infrastructure in Medell&iacute;n


*** =left

- Medellin, Colombia (population 2 million)
- As of early 2000s, 60% poverty rate, 20% unemployment, homicide 185 per 100K
- High residential segregation, w/ concentrated poverty in surrounding hills.
- 2004-2006: gondolas connect some but not all to city center.
- Cerda et al (2012, _Am J. Epideomiol._) study effects on neighborhood violence, propensity matching treatment neighborhoods to control.

*** =right

<embed height="625px" width="500px" src="./images/medellin-conc-pov.jpg"> 

*** =pnotes

- In contrast, Detroit's  2007 murder rate was 47/100K, and East St. Louis's was 102 per 100K.

--- &twocolleftwider


## Propensity scores in the Medellin study

*** =left

>- Small matching problem: $n_t=25$, $n_c=23$ (neighborhoods).
>- $\bar{x}_t$s and $\bar{x}_c$s not too far apart, but PS matching brought them closer.
>- Figure shows PS we matched on - the $X\hat{\beta}$ from a logistic regression (Cerda et al 2012, _Am J Epi_).
>- Region of strict overlap contains only 6/25 $t$s and 4/23 $c$s!
>- Yet Cerda et al matched all 48 neighborhoods.  (Using "full matching", which permits many-to-one matches.)

*** =right

<embed height="625px" width="400px" src="./images/boxplot-Medellin-logit.jpg">


--- 
# Outline

1. Problem: determining the region of  overlap
2. <span class = 'red'>Propensity score caliper matching</span>
3. PPSE calipers
4. Applications & discussion    





--- &twocolleftwider

## Matching with PS calipers (Rosenbaum & Rubin) 

>- For matching, useful to specify a tolerance or _caliper_ (Althauser & Rubin, 1970). 
>- Absent unmeasured confounding, matches with small PS differences mimic paired random assignment (R.&R. 1983). 
>- In effect, PS caliper (R.&R. 1985) adjudicates "closeness" to region of common support. 

*** =left

>- R. & R. (1985), Rubin & Thomas (2000) recommend $s_p(\widehat{\mathrm{PS}})/4$.
>- (Not $\widehat{\mathrm{PS}}$ the estimated probability, $\widehat{\mathrm{PS}}$ the estimated index function.
>  I.e., logits of estimated probabilities.)
>- The other widely used method (e.g. Austin & Lee 2009, Lunt 2013), besides requiring strict overlap.

*** =right

<embed height="400px" width="200px" src="./images/boxplot-Medellin-logit.jpg">



--- &twocol

## Room for improvement (in caliper = $.25s_p$)

>- $s_p$ approaches a positive limit, not 0. For large $n$, matches as crude as $.25s_p$ aren't very RCT-like.  (Even if there's no unmeasured confounding.) 
>- In the Medellin study (smallish $n$), $.25s_p(\widehat{\mathrm{PS}})$ caliper excludes 56% of treatment group.
>- At the same time, in that study it's tenable ($p=.10$) that all PSes are the same!

*** =left

>- <embed height="400px" width="400px" src="./images/boxplots-Medellin-logit.jpg">

*** =right

>- Moral: in small studies, $.25s_p(\widehat{\mathrm{PS}})$ may be too strict. 
>- In large studies, PS may be estimable w/ much more precision than $.25s_p(\widehat{\mathrm{PS}})$.  At same time, more potential controls. 
>- Moral: in large studies, $.25s_p(\widehat{\mathrm{PS}})$ may be too loose.

--- 

# Outline

1. Problem: determining the region of  overlap</span>
2. Propensity score caliper matching
3. <span class = 'red'>PPSE calipers</span>
4. Applications & discussion    


---

## Sampling variability of paired PS differences

>- Intuition: Even if we had matched **perfectly** on the PS, there would still be matched differences in $\widehat{\mathrm{PS}}$.  Let's estimate the size of these differences & use result to define caliper.    
>    - Should imply that as $n\uparrow \infty$, caliper $\downarrow 0$. Then, $\hat\beta \stackrel{P}{\rightarrow} \beta \Rightarrow$ should mean paired PS differences $\downarrow 0$.
>    - Excludes treatment subjects whose PSes are distinguishably different from all controls'.    
>- To limit model sensitivity, work on index scale.
>- If $\mathbf{d} = \mathbf{x}_1 - \mathbf{x}_2$, $\mathbf{d}(\hat{\beta} - \beta)$ is the error of estimation of paired PS difference $\mathbf{d}\beta$.
>- Considering all ${n \choose 2}$ possible pairs $r$, the expected MS estimation error of paired PS differences can be expressed as the Frobenius inner product, $\langle \cdot, \cdot \rangle_F$, of two covariances:   
\[ \mathbf{E}_\beta \left[{n \choose 2}^{-1}\sum_r (\mathbf{d}_r (\hat{\beta} - \beta))^2 \right] = 2\langle S^{(x)}, \mathbf{E}_\beta \{(\hat\beta - \beta) (\hat{\beta} -\beta)' \}\rangle_F . \]

---

## Sampling variability of paired PS diffs (ii)

>- But **all possible** pairs include many that are quite different on the PS.  We wanted to characterize $\mathbf{E} (\mathbf{d}_r(\hat{\beta} - \beta))^2$ among pairs $r$ s.t. $\mathbf{d}_r\beta \approx 0$.  
>- Since there may be few or no such pairs, instead "residualize" all pairs for the true PS, leaving remaining differences in place.  For $i=1,\ldots, p$, define $d_{(i)}^\perp =$ residual of $d_{(i)}$ regression on $\mathbf{d}\beta$; $\mathbf{d}^\perp = (d^\perp_1\,d^\perp_2\, \ldots\, d^\perp_p)$.  (By construction, $\mathbf{d}_r^\perp \beta =0$ for all $r$.)  

>- Rather than the mean-square of PS-difference estimation errors, ${n \choose 2}^{-1} \mathbf{E} \sum_r (\mathbf{d}_r(\hat{\beta} - \beta))^2$, consider the mean-square PS-difference estimation error, net of differences in the true PS: 
\[
{n \choose 2}^{-1}\mathbf{E}_\beta \sum_r (\mathbf{d}_r^\perp (\hat{\beta} - \beta))^2 =
2 \langle S^\perp, C_\hat\beta \rangle_F,
\]
where $S^\perp = \frac{1}{2}\mathrm{Cov}(\mathbf{d}^\perp) = \mathrm{Cov}(\mathbf{x}^{\perp \mathbf{x}\beta})$ and $x^{\perp \mathbf{x}\beta} = e(x |\mathbf{1}, \mathbf{x}\beta)$. 
>- Plugging in $\hat{S}^\perp = \mathrm{Cov}(\mathbf{x}^{\perp \mathbf{x}\hat{\beta}})$ and an estimate of ${C}_\hat\beta$ gives a squared "propensity-paired standard error" (PPSE).

---

# Outline

1. Problem: determining the region of  overlap
2. Propensity score caliper matching
3. PPSE calipers
4. <span class = 'red'>Applications & discussion</span>

---

## Example 2: Real or cosmetic reform in police practices?

Vagrancy arrests in the 60s and 70s

>- Before 1960, arrests for vagrancy were relatively common in U.S.; by 1980, relatively rare.
>- Increasingly forbidden by laws, court rulings, administrative decisions.
>- Did pretextual arrests decline, a civil rights victory, or were they simply shifted to other categories (displacement)?
>- Using administrative records assembled by R. Goluboff, D. Thacher and I identified 121 (local agency, year) instances, of $n=5600$ at-risk agency-years, that experienced a sharp, anomalaous decline in vagrancy arrests.
>- ... and we matched them to otherwise similar agencies in terms of year, region of U.S., population, prior use of vagrancy arrests and a propensity score. 


---

## Applications


>- For the Medellin example, 2.5 PPSEs works out to 4 logits ($3.9s_p$). 
>- No neighborhoods excluded. Recall that the alternatives excluded at least 12!
>- In the vagrancy study, 2.5 PPSEs = 2.3 logits ($1.9s_p$). 
>- Excludes a modest number of "treatment" cases (12/121).
>- Larger sample $\Rightarrow$ smaller PPSE. 
>- Even in larger samples I've worked with, tends to be more permissive than $.25s_p$ rule.

*** =pnotes

- In vagrancy study, PPSE criterion excluded 12/121 Txes, wheres .2$s_p$ would have excluded more.



---

## Calipers to ensure that PS differences tend to 0

>- Pairings (denoted "$i \sim j$") should satisfy $\max_{i \sim j} |(\vec{x}_i - \vec{x}_j)\beta| \downarrow 0$ as $n \uparrow \infty$.
>- Can this be accomplished with a requirement $\max_{i \sim j} |(\vec{x}_i - \vec{x}_j)\hat\beta| \leq w_n$, some $w_n \downarrow 0$?  (In light of $|(\vec{x}_i - \vec{x}_j)\beta| \leq |(\vec{x}_i - \vec{x}_j)(\hat\beta - \hat\beta)| + |(\vec{x}_i - \vec{x}_j)\hat\beta|$.)
>- ($w_n = s_p/4$ would mean $w_n \stackrel{P}{\rightarrow} \text{const} > 0$. Again, this won't do; we need $w_n \stackrel{P}{\rightarrow} 0$.)
>- For any $i,j\leq n$, $|(\vec{x}_i - \vec{x}_j)(\hat\beta - \beta)| \leq |\vec{x}_i - \vec{x}_j|_2|\hat\beta - \beta|_2$.
>- Expect $\max_{i,j}|(\vec{x}_i - \vec{x}_j)(\hat\beta - \beta)| \approx (\max_{i,j}|\vec{x}_i - \vec{x}_j|_2)|\hat\beta - \beta|_2$.
>- Even with bounded $x$s, this is $O(p^{1/2})O_P\left[(p/n)^{1/2}\right] = O_P(p/n^{1/2})$ (He & Shao, 2000).
>- We'll need to **assume** $p/n^{1/2} \downarrow 0$ --- not only $n \gg p$, also  $n \gg p^2$! 
>- A similar condition ($p^2 \log(p)/n \downarrow 0$) suffices for validity $\hat{C}_{\hat\beta}$, the Huber-White estimate of $\mathrm{Cov} (\hat\beta)$.  And for $\text{PPSE} = \left[2\langle \hat{S}^\perp, \hat{C}_{\hat\beta}\rangle_F\right]^{1/2}  \downarrow 0$ as $n \uparrow \infty$. 

---


### Discussion

>- Roughly speaking, the PPSE is the r.m.s. of $\vec{x}_i\hat\beta-\vec{x}_j\hat\beta$ among pairs $(i,j)$ s.t. $\vec{x}_i\beta = \vec{x}_j\beta$. 
>- If $p^2/n \downarrow 0$, then PPSE $\stackrel{P}{\rightarrow} 0$, and matching within $k$ PPSEs suffices for $\max_{i\sim j} |\vec{x}_i\beta-\vec{x}_j\beta| \stackrel{P}{\rightarrow} 0$.
>- If $p^2/n \not\rightarrow 0$, no caliper imposed on $\vec{x}_i\hat\beta-\vec{x}_j\hat\beta$ entails $\max_{i\sim j} |\vec{x}_i\beta-\vec{x}_j\beta| \stackrel{P}{\rightarrow} 0$.
>- If your matches are farther than a few PPSEs, you're matching outside region of overlap. 
>- I like calipers of $k=2.5$ PPSEs.


<!-- Matches between 1 and 3 PPSEs might be considered "marginal" (Austin & Lee, 2009).-->

### Recommendations: 

>- I'm using sandwich estimates of $\mathrm{Cov} (\hat\beta) $, w/ some special sauce (cf. github.com/benthestatistician/PISE) to limit numerical instability.
>- Our `optmatch` R package (Fredrickson et al, 2016) does optimal pair and full matching (Gu & Rosenbaum, 1993; Hansen & Klopfer, 2006; Stuart & Green, 2008), readily accommodating PS calipers.

