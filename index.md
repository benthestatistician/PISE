---
title       : Propensity score calipers and the overlap condition
author      : Ben Hansen, UMich Statistics (bbh@umich.edu).  Project site - 
date        : IMA Precision Medicine workshop, September 2017
job         : github.com/benthestatistician/PISE (code/manuscript), benthestatistician.github.io/PISE (slides)
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
>- Rosenbaum & Rubin's (1983) solution is match on $\widehat{\mathrm{Pr}}(Z=1| \mathbf{X}=\mathbf{x})$, the PS.
>- Remaining problem: In a cell $\mathbf{x}_0$ s.t. $\mathrm{Pr}(Z=1| \mathbf{X}=\mathbf{x}_0)=1$, there's no data to estimate $\mathrm{Pr}(Y_c \in \cdot |  \mathbf{X}=\mathbf{x}_0)$, even as $n \uparrow \infty$.  (Whether or not you match or otherwise use PSes.)
>- The second problem leads people to reduce analytic sample to region of common support on $\widehat{\mathrm{PS}}$ (the method of "strict overlap").

--- &twocol
##  Example 1: violence & public infrastructure in Medell&iacute;n


*** =left

- Medellin, Colombia (population 2 million)
- As of early 2000s, 60\% poverty rate, 20\% unemployment, homicide 185 per 100K
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
>- Yet Cerda et al full-matched all 48 neighborhoods.

*** =right

<embed height="625px" width="400px" src="./images/boxplot-Medellin-logit.jpg">


--- 
# Outline

1. Problem: determining the region of  overlap
2. <span class = 'red'>Propensity score caliper matching</span>
3. PPSE calipers
4. Applications & discussion    





--- &twocolleftwider

## Matching within PS calipers (R.& R.)

>- For matching, useful to specify a tolerance or _caliper_. 
>- Absent unmeasured confounding, matches with small PS differences mimic paired random assignment. 
>- In effect, caliper adjudicates "closeness" to region of common support. 

*** =left

>- R. & R. (1985, _Amer. Statist._), Rubin & Thomas (2000, _JASA_) recommend $s_p(\widehat{\mathrm{PS}})/4$.
>- (Not $\widehat{\mathrm{PS}}$ the estimated probability, $\widehat{\mathrm{PS}}$ the estimated index function.
>  I.e., logits of estimated probabilities.)
>- The other widely used method (e.g. Austin & Lee 2009, Lunt 2013), besides requiring strict overlap.

*** =right

<embed height="400px" width="200px" src="./images/boxplot-Medellin-logit.jpg">



--- &twocol

## Room for improvement (in caliper = $.25s_p$)


>- In the Medellin study, $.25s_p(\widehat{\mathrm{PS}})$ caliper still excludes 14/25 treatment neighborhoods.
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

## Calipers as insurance that PS differences tend to 0

>- Pairings (denoted "$i \sim j$") should satisfy $\sup_{i \sim j} |(\vec{x}_i - \vec{x}_j)\beta| \downarrow 0$ as $n \uparrow \infty$.
>- Can this be accomplished with a requirement of form $|(\vec{x}_i - \vec{x}_j)\hat\beta| \leq w_n$, since $|(\vec{x}_i - \vec{x}_j)\beta| \leq |(\vec{x}_i - \vec{x}_j)(\hat\beta - \hat\beta)| + |(\vec{x}_i - \vec{x}_j)\hat\beta|$?
>- ($w_n = s_p/4$ won't do.  We need $w_n \stackrel{P}{\rightarrow} 0$.  In itself, even this won't be quite enough.)
>- Scanning pairs $i,j\leq n$, $|(\vec{x}_i - \vec{x}_j)(\hat\beta - \beta)| \leq |\vec{x}_i - \vec{x}_j|_2|\hat\beta - \beta|_2$.
>- Expect $\sup_{i,j}|(\vec{x}_i - \vec{x}_j)(\hat\beta - \beta)| \approx (\sup_{i,j}|\vec{x}_i - \vec{x}_j|_2)|\hat\beta - \beta|_2$.
>- Even w/ $|x_{ij}|$ bounded, uniformly in $i,j$, and $n$, this is $O(p^{1/2})O_P([p/n]^{1/2}) = O_P(p/n^{1/2})$.  We'll need to assume $p/n^{1/2} \downarrow 0$. 
>- In that case, a $w_n$ that's $O_P(p/n^{1/2})$ will do the trick.
>- E.g., a $w_n$ measuring average size of $|(\vec{x}_i - \vec{x}_j)(\hat\beta - \beta)|$.

---

## Sampling variability of paired PS differences

>- Intuition: Even if we had matched **perfectly** on the PS, there would still be matched differences in $\widehat{\mathrm{PS}}$.  Let's estimate the size of these differences & use result to define caliper.    
>    - caliper $\downarrow 0$ as $n\uparrow \infty$. So, $\hat\beta {\rightarrow} \beta \Rightarrow$ paired PS differences $\downarrow 0$.
>    - excludes treatment subjects whose PSes are distinguishably different from all controls'.    
>- To limit model sensitivity, work on index scale.
>- If $\mathbf{d} = \mathbf{x}_1 - \mathbf{x}_2$, $\mathbf{d}(\hat{\beta} - \beta)$ is the error of estimation of paired PS difference $\mathbf{d}\beta$.
>- Considering all ${n \choose 2} $ possible pairs $r$, the expected MS estimation error of paired PS differences can be expressed as a Frobenius inner product:   
\[ \mathbf{E}_\beta \sum_r (\mathbf{d}_r (\hat{\beta} - \beta))^2/{n \choose 2} = 2\langle S^{(x)}, \mathbf{E}_\beta \{(\hat\beta - \beta) (\hat{\beta} -\beta)' \}\rangle_F . \]

---

## Sampling variability of paired PS diffs (ii)

>- But **all possible** pairs include many that are quite different on the PS.  We wanted to characterize $\mathbf{E} (\mathbf{d}_r(\hat{\beta} - \beta))^2$ among pairs $r$ s.t. $\mathbf{d}_r\beta \approx 0$.  
>- Since there may be few or no such pairs, instead "residualize" all pairs for the true PS, leaving remaining differences in place.  For $i=1,\ldots, p$, define $d_{(i)}^\perp =$ residual of $d_{(i)}$ regression on $\mathbf{d}\beta$; $\mathbf{d}^\perp = (d^\perp_1\,d^\perp_2\, \ldots\, d^\perp_p)$.  By construction, $\mathbf{d}_r^\perp \beta =0$ for all $r$.  

>- Rather than the mean-square of PS-difference estimation errors, ${n \choose 2}^{-1} \mathbf{E} \sum_r (\mathbf{d}_r(\hat{\beta} - \beta))^2$, consider the mean-square PS-difference estimation error, net of differences in the true PS: 
\[
{n \choose 2}^{-1}\mathbf{E}_\beta \sum_r (\mathbf{d}_r^\perp (\hat{\beta} - \beta))^2 =
2 \langle S^\perp, C_\hat\beta \rangle_F,
\]
where $S^\perp = \frac{1}{2}\mathrm{Cov}(\mathbf{d}^\perp) = \mathrm{Cov}(\mathbf{x}^{\perp \mathbf{x}\beta})$, $x^{\perp \mathbf{x}\beta} = e(x |\mathbf{1}, \mathbf{x}\beta)$. 
>- Plugging in $\hat{S}^\perp = \frac{1}{n-1}e(\mathbf{x} | 1, \mathbf{x}\hat{\beta})'e(\mathbf{x} | 1, \mathbf{x}\hat{\beta}) )$ and a $\hat{C}_\hat\beta$  extracted from the regression fit gives the "propensity-paired standard error" (PPSE).

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



--- &twocol

## Example 3: Vascular closure devices vs manual closure following percutaneous coronary intervention


*** =left

>- Once this stent has been threaded up into your heart, the hole in your femoral artery needs to be closed. 
>- VCDs are more comfortable than manual closure -- are they as safe?
>- Gurm et al (2013, _Ann Intern Med_) used data from a 32-hospital collaborative in Michgan to conduct a propensity-matched study. 
>- Large matching problem: $n_t=31$K; $n_c=54$K.  

*** =right

<embed height="400px" width="400px" src="./images/Stent-PCI.jpg">



---

## Applications


>- For the Medellin example, 2.5 PPSEs works out to 4 logits ($3.9s_p$). 
>- No neighborhoods excluded. Recall that the alternatives excluded at least 12!
>- In the vagrancy study, 2.5 PPSEs = 3.8 logits ($.5s_p$). 
>- In the VCD example 2.5 PPSEs= 0.15 logits ($.3s_p$). 
>- Larger sample $\Rightarrow$ smaller PPSE. 

*** =pnotes

- In vagrancy study, PPSE criterion excluded 2/121 Txes, wheres .2$s_p$ would have excluded 9.

--- &twocol

## Multiple PS calipers in the vagrancy study (Ex 2)


*** =left

<embed height="450px" width="450px" src="./images/boxplot-vagrancy-1ststage.jpg">

*** =right

<embed height="450px" width="450px" src="./images/boxplot-vagrancy-2ndstage.jpg">

---



### Discussion

>- Getting consistency out of PPSE caliper requires $p^2/n \downarrow 0$.  If $p/n \downarrow 0$, still interpretable as an SE of sorts. 
>- No need to restrict yourself to a single index score.
>- With a saturated PS model, observed information may be poorly conditioned.  The method of PPSE estimation as described here requires some elaboration (see code on github).
>- Likeness of observations ("like to like") is interpreted in terms of PS variables.  If too strict or too loose, then adjust PS variables.  

<!-- Matches between 1 and 3 PPSEs might be considered "marginal" (Austin & Lee, 2009).-->

### Recommendations: 

>- I like calipers of 2.5 PPSEs
>- I'm using sandwich estimates of $\mathrm{Cov} (\hat\beta) $, w/ some special sauce to limit numerical instability.
>- Our `optmatch` R package (Fredrickson et al, 2016) does optimal pair and full matching (Gu & Rosenbaum, 1993; Hansen & Klopfer, 2006; Stuart & Green, 2008), readily accommodating PS calipers.

