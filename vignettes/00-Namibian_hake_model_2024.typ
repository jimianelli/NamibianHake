// Some definitions presupposed by pandoc's typst output.
#let blockquote(body) = [
  #set text( size: 0.92em )
  #block(inset: (left: 1.5em, top: 0.2em, bottom: 0.2em))[#body]
]

#let horizontalrule = [
  #line(start: (25%,0%), end: (75%,0%))
]

#let endnote(num, contents) = [
  #stack(dir: ltr, spacing: 3pt, super[#num], contents)
]

#show terms: it => {
  it.children
    .map(child => [
      #strong[#child.term]
      #block(inset: (left: 1.5em, top: -0.4em))[#child.description]
      ])
    .join()
}

// Some quarto-specific definitions.

#show raw.where(block: true): block.with(
    fill: luma(230), 
    width: 100%, 
    inset: 8pt, 
    radius: 2pt
  )

#let block_with_new_content(old_block, new_content) = {
  let d = (:)
  let fields = old_block.fields()
  fields.remove("body")
  if fields.at("below", default: none) != none {
    // TODO: this is a hack because below is a "synthesized element"
    // according to the experts in the typst discord...
    fields.below = fields.below.amount
  }
  return block.with(..fields)(new_content)
}

#let empty(v) = {
  if type(v) == "string" {
    // two dollar signs here because we're technically inside
    // a Pandoc template :grimace:
    v.matches(regex("^\\s*$")).at(0, default: none) != none
  } else if type(v) == "content" {
    if v.at("text", default: none) != none {
      return empty(v.text)
    }
    for child in v.at("children", default: ()) {
      if not empty(child) {
        return false
      }
    }
    return true
  }

}

// Subfloats
// This is a technique that we adapted from https://github.com/tingerrr/subpar/
#let quartosubfloatcounter = counter("quartosubfloatcounter")

#let quarto_super(
  kind: str,
  caption: none,
  label: none,
  supplement: str,
  position: none,
  subrefnumbering: "1a",
  subcapnumbering: "(a)",
  body,
) = {
  context {
    let figcounter = counter(figure.where(kind: kind))
    let n-super = figcounter.get().first() + 1
    set figure.caption(position: position)
    [#figure(
      kind: kind,
      supplement: supplement,
      caption: caption,
      {
        show figure.where(kind: kind): set figure(numbering: _ => numbering(subrefnumbering, n-super, quartosubfloatcounter.get().first() + 1))
        show figure.where(kind: kind): set figure.caption(position: position)

        show figure: it => {
          let num = numbering(subcapnumbering, n-super, quartosubfloatcounter.get().first() + 1)
          show figure.caption: it => {
            num.slice(2) // I don't understand why the numbering contains output that it really shouldn't, but this fixes it shrug?
            [ ]
            it.body
          }

          quartosubfloatcounter.step()
          it
          counter(figure.where(kind: it.kind)).update(n => n - 1)
        }

        quartosubfloatcounter.update(0)
        body
      }
    )#label]
  }
}

// callout rendering
// this is a figure show rule because callouts are crossreferenceable
#show figure: it => {
  if type(it.kind) != "string" {
    return it
  }
  let kind_match = it.kind.matches(regex("^quarto-callout-(.*)")).at(0, default: none)
  if kind_match == none {
    return it
  }
  let kind = kind_match.captures.at(0, default: "other")
  kind = upper(kind.first()) + kind.slice(1)
  // now we pull apart the callout and reassemble it with the crossref name and counter

  // when we cleanup pandoc's emitted code to avoid spaces this will have to change
  let old_callout = it.body.children.at(1).body.children.at(1)
  let old_title_block = old_callout.body.children.at(0)
  let old_title = old_title_block.body.body.children.at(2)

  // TODO use custom separator if available
  let new_title = if empty(old_title) {
    [#kind #it.counter.display()]
  } else {
    [#kind #it.counter.display(): #old_title]
  }

  let new_title_block = block_with_new_content(
    old_title_block, 
    block_with_new_content(
      old_title_block.body, 
      old_title_block.body.body.children.at(0) +
      old_title_block.body.body.children.at(1) +
      new_title))

  block_with_new_content(old_callout,
    new_title_block +
    old_callout.body.children.at(1))
}

// 2023-10-09: #fa-icon("fa-info") is not working, so we'll eval "#fa-info()" instead
#let callout(body: [], title: "Callout", background_color: rgb("#dddddd"), icon: none, icon_color: black) = {
  block(
    breakable: false, 
    fill: background_color, 
    stroke: (paint: icon_color, thickness: 0.5pt, cap: "round"), 
    width: 100%, 
    radius: 2pt,
    block(
      inset: 1pt,
      width: 100%, 
      below: 0pt, 
      block(
        fill: background_color, 
        width: 100%, 
        inset: 8pt)[#text(icon_color, weight: 900)[#icon] #title]) +
      block(
        inset: 1pt, 
        width: 100%, 
        block(fill: white, width: 100%, inset: 8pt, body)))
}



#let article(
  title: none,
  authors: none,
  date: none,
  abstract: none,
  cols: 1,
  margin: (x: 1.25in, y: 1.25in),
  paper: "us-letter",
  lang: "en",
  region: "US",
  font: (),
  fontsize: 11pt,
  sectionnumbering: none,
  toc: false,
  toc_title: none,
  toc_depth: none,
  toc_indent: 1.5em,
  doc,
) = {
  set page(
    paper: paper,
    margin: margin,
    numbering: "1",
  )
  set par(justify: true)
  set text(lang: lang,
           region: region,
           font: font,
           size: fontsize)
  set heading(numbering: sectionnumbering)

  if title != none {
    align(center)[#block(inset: 2em)[
      #text(weight: "bold", size: 1.5em)[#title]
    ]]
  }

  if authors != none {
    let count = authors.len()
    let ncols = calc.min(count, 3)
    grid(
      columns: (1fr,) * ncols,
      row-gutter: 1.5em,
      ..authors.map(author =>
          align(center)[
            #author.name \
            #author.affiliation \
            #author.email
          ]
      )
    )
  }

  if date != none {
    align(center)[#block(inset: 1em)[
      #date
    ]]
  }

  if abstract != none {
    block(inset: 2em)[
    #text(weight: "semibold")[Abstract] #h(1em) #abstract
    ]
  }

  if toc {
    let title = if toc_title == none {
      auto
    } else {
      toc_title
    }
    block(above: 0em, below: 2em)[
    #outline(
      title: toc_title,
      depth: toc_depth,
      indent: toc_indent
    );
    ]
  }

  if cols == 1 {
    doc
  } else {
    columns(cols, doc)
  }
}

#set table(
  inset: 6pt,
  stroke: none
)
#show: doc => article(
  title: [Namibian hake model update, 2024],
  authors: (
    ( name: [John Kathena, Jim Ianelli],
      affiliation: [],
      email: [] ),
    ),
  date: [2024-07-18],
  toc: true,
  toc_title: [Contents],
  toc_depth: 3,
  cols: 1,
  doc,
)


#block[
```r
library(knitr)
hook_output <- knit_hooks$get("output")
knit_hooks$set(output = function(x, options) {
  # this hook is used only when the linewidth option is not NULL
  if (!is.null(n <- options$linewidth)) {
    x <- knitr:::split_lines(x)
    # any lines wider than n should be wrapped
    if (any(nchar(x) > n)) x <- strwrap(x, width = n)
    x <- paste(x, collapse = "\n")
  }
  hook_output(x, options)
})

knitr::opts_chunk$set(collapse = TRUE, comment = "  ", fig.align = "center", cache = FALSE, tidy.opts = list(width.cutoff = 80), tidy = TRUE)
knitr::opts_knit$set(root.dir = here::here())
# knitr::opts_chunk$set(warning=F, message=F, echo=F, results=F,fig.width=6, fig.height=5)
```

]
```yaml
project:
  render: ['*.qmd']
     - name: Set up Quarto
       uses: quarto-dev/quarto-actions/setup@v2
```

= Project overview
<project-overview>
The following reflects my interpretation of the Marine Stewardship Council’s certification request. There is a need to catch up on missed milestones and outlines the necessary steps for the upcoming Year 4 milestone.

=== Key Points
<key-points>
+ Year 3 Milestone Missed: The milestone required revised stock assessments for #emph[\M. paradoxus];, which was not met.
+ Year 4 Milestone: By February 2025, the MFMR must use the Harvest Control Rule (HCR) systematically to verify the TAC for #emph[\M. paradoxus];. This needs to be applied in the August/September 2024 management meetings.
+ Namibian Stock Assessment: There’s a recommendation to review and re-evaluate the assumptions and parameter values of assessment models, particularly the pessimistic base case model.
+ Implementation Issues: MFMR has Dr.~Ianelli’s report but not the code to run the model, and training is required for the Namibian team.

== Draft agenda
<draft-agenda>
The following draft agenda outlines the steps to address the missed milestones and prepare for the upcoming Year 4 milestone.

=== Week 1:
<week-1>
+ Day 1-2: Review and Planning
  - Review the Year 3 milestone requirements and current progress.
  - Plan steps to implement the HCR for #emph[\M. paradoxus];.
+ Day 3-4: Data Preparation
  - Gather and prepare Namibia stock assessment data.
  - Coordinate with MFMR to understand current data handling and management practices.
+ Day 5: Meeting Preparation
  - Prepare documentation and a presentation for the MFMR management meeting.
  - Outline the steps needed for the August/September 2024 meeting to include HCR in TAC setting.

=== Week 2:
<week-2>
#block[
#set enum(numbering: "1.", start: 4)
+ Day 1-2: Model Review
  - Review developments and report key elements for implementation.
  - Develop a preliminary implementation plan for the HCR model.
+ Day 3-4: Training Coordination
  - Arrange a training session with Dr.~Ianelli or another suitable individual for MFMR.
  - Coordinate with the training provider and MFMR to schedule the session.
+ Day 5: Reporting
  - Compile a progress report summarizing activities, challenges, and next steps.
  - Send the report to Hugh and relevant stakeholders for feedback.
]

This agenda ensures a systematic approach to address the milestones and prepare for the upcoming management meeting, focusing on implementing the HCR and providing necessary training to MFMR.

Below are two main sections, first on model developments and second on application of the control rule that accounts for the signals in the data on the different species.

= Assessment model runs
<assessment-model-runs>
The original base-case model was evaluated for a number of features and extensions. These included focus on what data components were fit well and how improvements in consistency can be made. For the latter part, we found that the fits to the index and CPUE data were particularly poor and could be improved. We reviewed the model from Ianelli et al.~(2023) and decided that whilst reassuring that alternative assessment approach provide similar results, for consistency and due to the added features of the Namibian model, it was preferred to make necessary changes in that code base. Give the short time available for training and implementation, we decided to focus on the base-case model and make necessary changes. The following sections outline the model runs and the results relative to the previous assessment.

=== Model descriptions
<model-descriptions>
The following table was developed based on testing the model with different assumptions and data sources. Key differences from the 2023 assessment configuration was the assumption that model estimation of variance terms was appropriate. This feature resulted in unacceptable residual patterns and essentially a complete down weighting of the index data. We used the assumed variance terms (CVs) for the indices in all of the following model configurations:

#figure(
  align(center)[#table(
    columns: (33.8%, 66.2%),
    align: (auto,auto,),
    table.header([Model], [Description],),
    table.hline(),
    [Previous base case], [As specified in past assessments, estimated steepness and all variance terms],
    [Base case (m0)], [Model with survey "minus group" to be ages 0, and 1 instead of 0, 1, and 2 as done in the past, steepness fixed at 0.7, q estimated, and time-varying fishery asymptotic selectivity specified.],
    [m1], [As base case but with survey catchability fixed at 1.0],
    [m2], [As base case but with survey catchability fixed at 0.5],
    [m3], [As base case but with natural mortality estimated],
    [m4], [As base case but with fishery selectivity allowed to be dome-shaped],
    [m5], [As base case but with stock-recruit steepness fixed at 0.5],
    [m6], [As base case but with stock-recruit steepness fixed at 0.9],
  )]
  , kind: table
  )

#block[
```r
library(NamibianHake)
library(flextable)
library(here)
library(tidyverse)
library(ggridges)
theme_set(ggthemes::theme_few())
library(xtable)
library(kableExtra)
set_flextable_defaults(digits = 3, decimal.mark = ".", big.mark = ",", na_str = "<na>")
```

]
#block[
```r
mod_ref <- c("old_bc", "m0", "m1", "m2", "m3", "m4", "m5", "m6")
mod_dir <- c("old_bc", "m0", "m1", "m2", "m3", "m4", "m5", "m6")
mod_label <- c("2023 base case", "2024 base case", "Model 1", "Model 2", "Model 3",
    "Model 4", "Model 5", "Model 6")


#---Main code that extracts all the results from the model lists above------
res <- get_results(mod_names. = mod_label, moddir = mod_dir)

modlst <- res$modlst
old_bc <- modlst[[1]]
m0 <- modlst[[2]]
moddiag <- res$moddiag
dfsrr <- data.frame()
for (i in 1:length(mod_ref)) {
    dfsrr <- rbind(dfsrr, data.frame(Model = mod_label[i], SSB = modlst[[i]]$SSB,
        R = modlst[[i]]$Pred_Rec))
}
mods <- data.frame()
for (i in 1:length(mod_ref)) {
    mods <- rbind(mods, data.frame(moddiag[[i]], Model = names(moddiag[i])))
}
```

]
=== Fits to index data
<fits-to-index-data>
After delving into the details of the model specifications, our main conclusion was that the previous assessments were generally insensitive to the trend data (surveys and CPUE series). This was due to the fact that the variance terms (CVs) were estimated. A more standard approach, where CVs are specified based on the data (e.g., from design-based sampling theory), was used and this performed better (@fig-fitidx1, @fig-fitcpue3). In the previous assessment, the fit tot he early period of CPUE data was better, but in this case the CV was estimated to be about 10%, and extremely low value for this type of data (@fig-fitcpue1).

Similarly, we re-evaluated the ability of this model to estimate the stock recruitment productivity parameter (steepness). It is extremely rare that sufficient data are available to freely estimate this, even with extensive high-quality index data. Therefore we evaluated what assumptions were taken elsewhere (for this and related species) and ran the models with steepness fixed at 0.7.

#figure([
#box(image("00-Namibian_hake_model_2024_files/figure-typst/fig-fitidx1-1.svg"))
], caption: figure.caption(
position: bottom, 
[
Base case model fits to main survey data compared to the previous assessment.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-fitidx1>


#figure([
#box(image("00-Namibian_hake_model_2024_files/figure-typst/fig-fitcpue3-1.svg"))
], caption: figure.caption(
position: bottom, 
[
Base case model fits to the CPUE index 3 data compared to the previous assessment.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-fitcpue3>


```r
dfcpue <- rbind(data.frame(Year = 1964:2023, Obs = old_bc$Obs_CPUE_1, predicted = old_bc$e_CPUE_1,
    Model = "Base case 2023"), data.frame(Year = 1964:2024, Obs = m0$Obs_CPUE_1,
    predicted = m0$e_CPUE_1, Model = "Base case 2024"))
dfcpue |>
    filter(Obs > 0) |>
    ggplot(aes(x = Year, y = Obs, color = Model)) + geom_point(color = "black") +
    geom_line(aes(y = predicted)) + geom_point(aes(y = predicted)) + ggtitle("CPUE 1") +
    ylim(0, NA)
ggsave(here("mods", "figs", "cpue1.png"), width = 6, height = 5)
```

#figure([
#box(image("00-Namibian_hake_model_2024_files/figure-typst/fig-fitcpue1-1.svg"))
], caption: figure.caption(
position: bottom, 
[
Base case model fits to the CPUE index 1 data compared to the previous assessment.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-fitcpue1>


=== SSB results
<ssb-results>
@fig-ssb shows the SSB estimates for the base case model compared to the previous assessment. The horizontal lines correspond to the Bmsy values for the separate models.

```r
dftmp <- data.frame(Model = mod_label[1], Bmsy = (modlst[[1]]$Bmsy))
dftmp <- rbind(dftmp, data.frame(Model = mod_label[2], Bmsy = (modlst[[2]]$Bmsy)))
mods |>
    filter(Model %in% mod_label[1:2], Year > 1960, Variable == "SSB") |>
    ggplot(aes(x = Year, y = value, ymin = ymin, ymax = ymax, type = Model, fill = Model)) +
    geom_ribbon(alpha = 0.24) + ggthemes::theme_few() + coord_cartesian(ylim = c(0,
    6000)) + geom_line(aes(color = Model)) + geom_hline(yintercept = dftmp$Bmsy[1],
    color = 2) + geom_hline(yintercept = dftmp$Bmsy[2], color = 3) + geom_point(aes(color = Model,
    shape = Model), size = 1) + ylab("SSB") + xlab("Year")
ggsave(here("mods", "figs", "ssb1.png"), width = 6, height = 5)
```

#figure([
#box(image("00-Namibian_hake_model_2024_files/figure-typst/fig-ssb-1.svg"))
], caption: figure.caption(
position: bottom, 
[
Base case model showing the SSB estimates compared to the previous assessment. The horizontal lines correspond to the Bmsy values for the separate models.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-ssb>


Another way to evaluate the relative trends in spawning biomass is to examine the so-called "depletion" levels. This is the ratio of the current SSB to the theoretical unfished value. The results are shown in @fig-depl for the base case model and in @fig-deplalt for the alternative models.

```r
mods |>
    filter(Model %in% mod_label[1:2], Year > 1960, Variable == "Depletion") |>
    ggplot(aes(x = Year, y = value, ymin = ymin, ymax = ymax, type = Model, fill = Model)) +
    geom_ribbon(alpha = 0.24) + ggthemes::theme_few() + geom_line(aes(color = Model)) +
    geom_point(aes(color = Model, shape = Model), size = 1) + ylab("Relative SSB") +
    xlab("Year")
ggsave(here("mods", "figs", "depl1.png"), width = 6, height = 5)
```

#figure([
#box(image("00-Namibian_hake_model_2024_files/figure-typst/fig-depl-1.svg"))
], caption: figure.caption(
position: bottom, 
[
Base case model showing the relative SSB estimates compared to the previous assessment.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-depl>


```r
mods |>
    filter(Model %in% mod_label[2:9], Year > 1960, Variable == "Depletion") |>
    ggplot(aes(x = Year, y = value, ymin = ymin, ymax = ymax, type = Model, fill = Model)) +
    geom_ribbon(alpha = 0.24) + ggthemes::theme_few() + geom_line(aes(color = Model)) +
    geom_point(aes(color = Model, shape = Model), size = 1) + ylab("Relative SSB") +
    xlab("Year")
   Warning: The shape palette can deal with a maximum of 6 discrete values because more
   than 6 becomes difficult to discriminate
   ℹ you have requested 7 values. Consider specifying shapes manually if you need
     that many have them.
   Warning: Removed 61 rows containing missing values or values outside the scale range
   (`geom_point()`).
ggsave(here("mods", "figs", "depl2.png"), width = 6, height = 5)
   Warning: The shape palette can deal with a maximum of 6 discrete values because more
   than 6 becomes difficult to discriminate
   ℹ you have requested 7 values. Consider specifying shapes manually if you need
     that many have them.
   Removed 61 rows containing missing values or values outside the scale range
   (`geom_point()`).
```

#figure([
#box(image("00-Namibian_hake_model_2024_files/figure-typst/fig-deplalt-1.svg"))
], caption: figure.caption(
position: bottom, 
[
Alternative model results showing the relative SSB estimates.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-deplalt>


The recruitment results are consistent with those seen for the spawning biomass. The base case model shows a decline in recruitment estimates compared to the previous assessment, as shown in @fig-rec. Compared to the base case model, the alternative models show a range of recruitment estimates, as shown in @fig-recalt.

```r
mods |>
    filter(Model %in% mod_label[1:2], Year > 2010, Variable == "R") |>
    ggplot(aes(x = Year, y = value, ymin = ymin, ymax = ymax, type = Model, fill = Model)) +
    geom_errorbar(width = 0.95, position = "dodge", alpha = 0.3) + geom_bar(width = 0.95,
    stat = "Identity", position = "dodge") + ggthemes::theme_few() + ylab("Recruitment age 0") +
    xlab("Year")
ggsave(here("mods", "figs", "rec.png"), width = 6, height = 5)
```

#figure([
#box(image("00-Namibian_hake_model_2024_files/figure-typst/fig-rec-1.svg"))
], caption: figure.caption(
position: bottom, 
[
Alternative model results on recruitment estimates.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-rec>


```r
mods |>
    filter(Model %in% mod_label[2:9], Year > 2010, Variable == "R") |>
    ggplot(aes(x = Year, y = value, ymin = ymin, ymax = ymax, type = Model, fill = Model)) +
    geom_errorbar(width = 0.95, position = "dodge", alpha = 0.3) + geom_bar(width = 0.95,
    stat = "Identity", position = "dodge") + ggthemes::theme_few() + ylab("Recruitment age 0") +
    xlab("Year")
ggsave(here("mods", "figs", "recalt.png"), width = 6, height = 5)
```

#figure([
#box(image("00-Namibian_hake_model_2024_files/figure-typst/fig-recalt-1.svg"))
], caption: figure.caption(
position: bottom, 
[
Alternative model results on recruitment estimates.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-recalt>


To show the history relative to the replacement yield, we can plot the SSB/Bmsy and the catch/replacement yield (@fig-kobe). This shows a significant difference between the previous assessment and that proposed as the base-case for this year. The alternative models show a range of results, as shown in @fig-kobe2.

```r
tmp <- mods |>
    select(Year, Model, Variable, value) |>
    filter(Model %in% mod_label[1:2], Year > 1980, Variable %in% c("Catch_RY", "B_Bmsy")) |>
    pivot_wider(names_from = Variable, values_from = value) |>
    arrange(Model, Year)  #|>
tmp |>
    ggplot(aes(x = B_Bmsy, label = Year, y = Catch_RY, shape = Model, color = Model,
        fill = Model)) + geom_path() + ggthemes::theme_few() + geom_point() + geom_text(alpha = 0.5) +
    geom_hline(yintercept = 1) + geom_vline(xintercept = 1) + xlab("SSB/Bmsy") +
    ylab("Catch / replacement Yield")
ggsave(here("mods", "figs", "kobe1.png"), width = 6, height = 5)
```

#figure([
#box(image("00-Namibian_hake_model_2024_files/figure-typst/fig-kobe-1.svg"))
], caption: figure.caption(
position: bottom, 
[
Base case model showing the relative SSB estimates compared to the previous assessment.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-kobe>


```r
tmp <- mods |>
    select(Year, Model, Variable, value) |>
    filter(Model %in% mod_label[1:2], Year > 1980, Variable %in% c("Catch_RY", "Depletion")) |>
    pivot_wider(names_from = Variable, values_from = value)
tmp <- mods |>
    select(Year, Model, Variable, value) |>
    filter(Model %in% mod_label[c(2:3, 7:8)], Year > 1980, Variable %in% c("Catch_RY",
        "B_Bmsy")) |>
    pivot_wider(names_from = Variable, values_from = value) |>
    arrange(Model, Year)  #|>
tmp |>
    ggplot(aes(x = B_Bmsy, label = Year, y = Catch_RY, shape = Model, color = Model,
        fill = Model)) + geom_path() + ggthemes::theme_few() + geom_point() + geom_text(alpha = 0.5) +
    geom_hline(yintercept = 1) + geom_vline(xintercept = 1) + xlab("SSB/Bmsy") +
    ylab("Catch / replacement Yield")
ggsave(here("mods", "figs", "kobe2.png"), width = 6, height = 5)
```

#figure([
#box(image("00-Namibian_hake_model_2024_files/figure-typst/fig-kobe2-1.svg"))
], caption: figure.caption(
position: bottom, 
[
Base case model showing the relative SSB estimates compared to the model alternatives assessment.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-kobe2>


=== Age composition fits
<age-composition-fits>
The age composition fits for the base case model are shown inf @fig-agecomps1 and @fig-agecomps2. The base case model uses a 'minus group' equal to '1' for the survey data and for the fishery it was set to 2 (as was done previously).

```r
PlotAgeFit(x = m0, title = "Base case", type = "fishery", fage = 2, lage = 7)
ggsave(here("mods", "figs", "Age_comp_fish.png"), width = 9, height = 8)
```

#figure([
#box(image("00-Namibian_hake_model_2024_files/figure-typst/fig-agecomps1-1.svg"))
], caption: figure.caption(
position: bottom, 
[
Base case model fits to the fishery age composition data. Note that the base case model uses a 'minus group' equal to '2' for the fishery data.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-agecomps1>


```r
PlotAgeFit(x = m0, title = "Base case", type = "survey1", fage = 1, lage = 7)
ggsave(here("mods", "figs", "Age_comp_surv.png"), width = 9, height = 8)
```

#figure([
#box(image("00-Namibian_hake_model_2024_files/figure-typst/fig-agecomps2-1.svg"))
], caption: figure.caption(
position: bottom, 
[
Base case model fits to the survey age composition data. Note that the base case model uses a 'minus group' equal to '1' for the survey data.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-agecomps2>


=== Selectivity
<selectivity>
#figure([
#box(image("00-Namibian_hake_model_2024_files/figure-typst/selex-1.svg"))
], caption: figure.caption(
position: bottom, 
[
Selectivity estimates for the base-case model run.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)


=== Stock-recruitment curves
<stock-recruitment-curves>
The following plot shows the stock-recruitment curves for the base case and models 1, 2, 3, and 5 and 6.

```r
p1 <- dfsrr |>
    ggplot(aes(x = SSB, y = R, color = Model, shape = Model)) + geom_point() + geom_line() +
    xlim(c(0, 4000))
p1
   Warning: The shape palette can deal with a maximum of 6 discrete values because more
   than 6 becomes difficult to discriminate
   ℹ you have requested 8 values. Consider specifying shapes manually if you need
     that many have them.
   Warning: Removed 47 rows containing missing values or values outside the scale range
   (`geom_point()`).
   Warning: Removed 7 rows containing missing values or values outside the scale range
   (`geom_line()`).
```

#figure([
#box(image("00-Namibian_hake_model_2024_files/figure-typst/srrplots-1.svg"))
], caption: figure.caption(
position: bottom, 
[
Model runs
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)


```r
ggsave(here("mods", "figs", "srr_curves.png"), width = 6, height = 5)
   Warning: The shape palette can deal with a maximum of 6 discrete values because more
   than 6 becomes difficult to discriminate
   ℹ you have requested 8 values. Consider specifying shapes manually if you need
     that many have them.
   Warning: Removed 47 rows containing missing values or values outside the scale range
   (`geom_point()`).
   Warning: Removed 7 rows containing missing values or values outside the scale range
   (`geom_line()`).
```

=== Stock status comparisons
<stock-status-comparisons>
The #cite(<tab-quant>, form: "prose") show the results of the base case model compared to the alternative models for a number of key statistics. Among these, the best fitting models were the base-case and Model 5. However, model 5 fits indicated that the

#block[
```r
# |
dftmp <- NULL
mod_scen <- c(2:8)
filler <- " "
names(filler) <- "-ln(Likelihood)"
for (ii in mod_scen) {
    x <- modlst[[ii]]
    nll <- round(x$ObjFun, 0)
    names(nll) <- "Overall"
    CPUE <- round(x$CPUE_Like, 0)
    names(CPUE) <- "CPUE"
    surv <- round(x$Survey_Like, 0)
    names(surv) <- "Survey"
    caa <- round(x$CAA_Likelihood, 0)
    names(caa) <- "Commercial CAA"
    caas <- round(x$CAAS_Likelihood, 0)
    names(caas) <- "Survey CAA"
    oneyrold <- round(x$Oneyearold_Likelihood, 0)
    names(oneyrold) <- "One yr-old biomass"
    rec <- round(x$RecRes_Likelihood, 0)
    names(rec) <- "Rec. resids."
    datalike <- nll - rec
    names(datalike) <- "Data likelihood sub-total"
    NumPars <- round(x$Npars, 0)
    names(NumPars) <- "Number parameters"
    AIC <- round(x$Akaike, 1)
    names(AIC) <- "AIC"
    v <- c(filler, nll, CPUE, surv, caa, caas, oneyrold, datalike, rec, NumPars,
        AIC)
    dftmp <- cbind(dftmp, v)
}
dftmp <- data.frame(rownames(dftmp), dftmp, row.names = NULL)
names(dftmp) <- c("Component", mod_label[mod_scen])
tabcap <- "Model fits to data components"
flextable(dftmp) |>
    set_caption(caption = tabcap) |>
    colformat_double() |>
    autofit()
```

#align(center)[
#box(image("00-Namibian_hake_model_2024_files/figure-typst/tab-quant-1.png"))
]
]
```r
dftmp <- NULL
mod_scen <- c(2:8)
for (ii in mod_scen) {
    x <- modlst[[ii]]
    Ksp <- round(x$KspSTD, 0)
    names(Ksp) <- "Unfished spawning biomass"
    Kexp <- round(x$KexpSTD, 0)
    names(Kexp) <- "Unfished expl. biomass"
    steepness <- round(x$Steep, 3)
    names(steepness) <- "SRR steepness"
    Cur_B <- round(x$Bstd[length(x$Bstd)], 3)
    names(Cur_B) <- "2024 SSB"
    Bmsy <- round(x$Spmsy, 0)
    names(Bmsy) <- "SSB_msy"
    Cur_B0 <- round(x$Cur_B0, 3)
    names(Cur_B0) <- "Current SSB over unfished"
    Cur_Bmsy <- round(x$Cur_Bmsy, 3)
    names(Cur_Bmsy) <- paste0("Current SSB over Bmsy")
    MSY <- round(x$MSY, 0)
    names(MSY) <- "MSY"
    ry <- round(x$aveRY_last5, 0)
    names(ry) <- paste0("recent 5-yr average replacement yield")
    ry_cur <- round(x$aveRY_last5 * x$Cur_Bmsy, 3)
    names(ry_cur) <- "Recent RY x current SSB /Bmsy"
    v <- c(Ksp, Kexp, steepness, Cur_B, Bmsy, Cur_B0, Cur_Bmsy, MSY, ry, ry_cur)
    dftmp <- cbind(dftmp, v)
}
dftmp <- data.frame(rownames(dftmp), dftmp, row.names = NULL)
names(dftmp) <- c("Statistic", mod_label[mod_scen])
tabcap <- "Selected management measures from alternative models. "
flextable(dftmp) |>
    set_caption(caption = tabcap) |>
    colformat_double() |>
    autofit()
```

#align(center)[
#box(image("00-Namibian_hake_model_2024_files/figure-typst/tab-stockstatus-1.png"))
]
```r
# align=paste0('lll',strrep('r',length(mod_scen+1)))) kable(tab,
# caption.placement = 'top', include.rownames = FALSE, sanitize.text.function =
# function(x){x}) print(tab)
```

A comparison among the models for stock status and

#figure([
#box(image("00-Namibian_hake_model_2024_files/figure-typst/fig-stockstatus-1.svg"))
], caption: figure.caption(
position: bottom, 
[
Some reference point comparisons among model runs. aveRY\_90 is the average replacement yeild since 1990, aveRY is the average replacement yield over the most recent 5 years before the current year, Cur\_90 is the current (terminal year) SSB over the estimate from 1990, Cur\_B0 is over the unfished estimate, and Cur\_Bmsy is the ratio of current SSB over the Bmsy estimate.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-stockstatus>


= Control rule application
<control-rule-application>
== Modeling Namibian hake survey data by species
<modeling-namibian-hake-survey-data-by-species>
Fisheries stock assessments require data that are reliably collected and compiled. Secondarily, assessment models should be configured to match the assumptions associated with the observed data. To account for survey trends between the two species of hake we applied the estimated observation errors to a simple state-space random walk model. This approach has a number of options for how process-errors can be specified and estimated. The observation model applies the observation-error variances $(sigma_(j , t)^2)$ for the $j^(t h)$ species in year $t (x_(j , t))$. The indices are fit to latent state variables, e.g., the underlying population trend $l n (hat(Z)_(j , t))$ as follows:

$l n (Z_(j , t)) = l n (hat(Z)_(j , t)) + epsilon.alt_(j , t)$ where $epsilon.alt_(j , t) tilde.op N (0 , sigma_(j , t)^2)$

and the state equation and associated process error variance $sigma_(P E)^2$ is defined as

$l n \( hat(Z)_(j , t + 1) = l n (hat(Z)_(j , t +)) + eta_(j , t) ,$ where $eta_(j , t) tilde.op N (0 , tau_j^2) .$

The process error variances $tau_j^2$ (which may or may not vary across indices) are fixed effect parameters and the unobserved species combined population $l n (Z_(j , t))$ is estimated as a series of random effects. The model is fit using maximum likelihood estimation in TMB using the R package "rema" (Sullivan 2022). The survey data for each species was used with CVs applied for observation error specifications. The values for $tau_j^2$ were tested for each species and found to be similar so they were set to the same values.

The above analysis provides a summary of the model runs and the design of a control rule that accounts for the signals in the data on the different species.

== Application to the Namibian hake stocks
<application-to-the-namibian-hake-stocks>
The control rule was first applied to the Namibian hake stocks using the survey data and the relative proportions of the two species. However, we wish to have a more robust control rule that accounts for the signals in the data on #emph[\M. paradoxus] biomass and have that be independent of #emph[\M. capensis] (e.g., @fig-relmean). For that case, we applied the survey smoothing model to #emph[paradoxus] alone. The next step was to compute mean biomass over the period 1990-2024, and evaluate the adustment for different levels of $gamma$. The historical adjustments based on this aspect of the MP is shown in @fig-adjust.

```
   Rows: 34 Columns: 4
   ── Column specification ────────────────────────────────────────────────────────
   Delimiter: ","
   chr (1): strata
   dbl (3): year, biomass, cv
   
   ℹ Use `spec()` to retrieve the full column specification for this data.
   ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
   Joining with `by = join_by(year, strata)`
   Model runtime: 0.1 seconds 
   stats::nlminb thinks the model has converged: mod$opt$convergence == 0
   Maximum gradient component: 6.86e-08 
   Max gradient parameter: log_PE 
   TMB:sdreport() was performed successfully for this model
   Joining with `by = join_by(strata, year)`
   Joining with `by = join_by(year)`
   Warning: Removed 2 rows containing missing values or values outside the scale range
   (`geom_point()`).
   Removed 2 rows containing missing values or values outside the scale range
   (`geom_point()`).
```

#figure([
#box(image("00-Namibian_hake_model_2024_files/figure-typst/fig-relmean-1.svg"))
], caption: figure.caption(
position: bottom, 
[
Historical biomass estimates for M. paradoxus and the mean biomass over the period 1990-2024.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-relmean>


#figure([
#box(image("00-Namibian_hake_model_2024_files/figure-typst/fig-adjust-1.svg"))
], caption: figure.caption(
position: bottom, 
[
Historical adjustments for different reactivities to #emph[\M. paradoxus] biomass.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-adjust>


== Control rule developments
<control-rule-developments>
The situation for developing a two-species control rule where catches between species and the trend in overall biomass for both species is combined is challenging. For management purposes, the goal is to avoid incidental takes in the proportion of one species that exceeds the historical levels of depletion for either stock. Fortunately, survey data are available that can be used to distinguish trends in the relative biomass for both stocks. The design of the triggered control rule therefore must consider patterns in the relative biomass from the survey data, the absolute biomass of the combined catch and biomass as modeled from the combined-stocks assessment. This provides a pragmatic approach using available data.

The steps in the control rule would be first to run a simple model that projected the survey biomass and relative proportion of #emph[\M. paradoxus];. Then, given the mean proportion over the period, compute the adjustment needed to the overarching control rule for the management procedure. For example, the historical range based on the survey has been between about one third of the mean value (in the earliest part of the period) to about 70% above the mean proportion. This range (especially the lower value) was used as a semi-empirical way to develop a minimum stock size threshold as part of the control rule. That is, when the stock of #emph[\M. paradoxus] drops below 30% of the mean proportion of the combined stocks estimate, the TAC recommendation for the combined stock would be zero. So if proportion of #emph[\M. paradoxus] $(p_y)$ is greater than 30% of the long-term mean proportion, then

$T A C_y = (R sum_y^(y - 4) frac(R Y_y, 5)) min (1.0 , B_y \/ B_(M S Y))^lambda 0.5^(- lambda) \) min (1 , dot(p_y) \/ dot(p^(‾)))^gamma$

where the rebuilding factor ($R$) is set to 0.8 when the spawning biomass is below $B_(M S Y)$ and 1.0 when above. In words, the TAC in year y is equal to the catch under the current control rule times the ratio of the spawning biomass relative to BMSY (or proxy) and the externally estimated proportion of #emph[\M. paradoxus] from survey data. The second term relates $B_(M S Y)$ and is intended to take fast action (for $lambda > 1.0$) when the biomass falls below 0.5 of that value (a standard in many places to define "overfished"). The third term on the right hand side reflects the impact the #emph[\M. paradoxus] biomass projected from a survey smoother (described below) and adjusts the TAC advice downwards when the projected $dot(p_y)$ drops below the mean value. The values of $gamma , lambda$ were evaluated and are shown in Table xx. We note that the specification of a survey linkage by the individual species provides an appropriate adjustment that reduces the exploitation rate and prevents potential for "the point of recruitment impairment" (PRI).

For the control rule as specified, the reactivity of the TAC advice to changes in either B/Bmsy or the biomass of #emph[\M. paradoxus] biomass can be adjusted. These are shown in @fig-gamma. For the purposes of this analysis, we set $gamma = .25$. For $lambda$, most models evaluated were above 0.5 of $B_(M S Y)$ so there was no added adjustment beyond the $R = 0.8$. So, following the TAC as specified, the base-case model was below \$B\_{MSY} so $R = 0.8$ with the average replacement yield over the last 5 years is recommended, with the adjustment based on the relative biomass of #emph[\M. paradoxus] and the long-term mean biomass of the species.

The TAC advice for the combined stocks is then given in the following table:

```r
dftmp <- NULL
mod_scen <- c(2:8)
ii = 2
names(cur_mp) <- "Ratio of paradoxus to mean "
for (ii in mod_scen) {
    x <- modlst[[ii]]
    Cur_Bmsy <- round(x$Cur_Bmsy, 3)
    names(Cur_Bmsy) <- paste0("Current B over Bmsy")
    ry <- round(x$aveRY_last5, 0)
    names(ry) <- paste0("recent 5-yr avg repl. yield")

    ry_cur <- round(x$aveRY_last5 * min(1, x$Cur_Bmsy), 3)
    # names(ry_cur) <- ('Recent RY x current B/Bmsy')

    lambda = 1.5
    gamma = 0.25
    if (Cur_Bmsy > 0.5)
        adj = 1 else adj = Cur_Bmsy^lambda/0.5^lambda
    tac1 <- round(0.8 * ry * adj * min(1, cur_mp)^gamma, 2)
    names(tac1) <- paste0("opt 1, lambda=", lambda, " gamma=", gamma)

    lambda = 2
    gamma = 0.25
    if (Cur_Bmsy > 0.5)
        adj = 1 else adj = Cur_Bmsy^lambda/0.5^lambda
    tac2 <- round(0.8 * ry * adj * min(1, cur_mp)^gamma, 2)
    names(tac2) <- paste0("opt 2, lambda=", lambda, " gamma=", gamma)

    lambda = 2
    gamma = 0.1
    if (Cur_Bmsy > 0.5)
        adj = 1 else adj = Cur_Bmsy^lambda/0.5^lambda
    adj
    Cur_Bmsy
    tac3 <- round(0.8 * ry * adj * min(1, cur_mp)^gamma, 2)
    names(tac3) <- paste0("opt 3, lambda=", lambda, " gamma=", gamma)

    lambda = 2
    gamma = 1
    if (Cur_Bmsy > 0.5)
        adj = 1 else adj = Cur_Bmsy^lambda/0.5^lambda
    adj
    Cur_Bmsy
    tac4 <- round(0.8 * ry * adj * min(1, cur_mp)^gamma, 2)
    names(tac4) <- paste0("opt 4, lambda=", lambda, " gamma=", gamma)

    filler = ""
    names(filler) <- "TAC"

    v <- c(ry, Cur_Bmsy, cur_mp, filler, tac1, tac2, tac3, tac4)
    dftmp <- cbind(dftmp, v)
}
dftmp <- data.frame(rownames(dftmp), dftmp, row.names = NULL)
names(dftmp) <- c("Statistic", mod_label[mod_scen])
tabcap <- "Selected management measures from alternative models. "
flextable(dftmp) |>
    set_caption(caption = tabcap) |>
    colformat_double() |>
    autofit()
```

#align(center)[
#box(image("00-Namibian_hake_model_2024_files/figure-typst/label-table-controlrule-1.png"))
]
````r
#--- Plot the adjustment factors for different reactivities to M. paradoxus biomass----
df01 <- tibble(relbiom = 1:20/20, `Gamma = 0.1` = (1:20/20)^0.1, `Gamma = 0.5` = (1:20/20)^0.5,
    `Gamma = 0.25` = (1:20/20)^0.25, ) |>
    pivot_longer(cols = -relbiom, names_to = "Reactivity", values_to = "adjustment")
df01 |>
    ggplot(aes(y = adjustment, x = relbiom, color = Reactivity)) + geom_line() +
    xlab("Relative M. paradoxus biomass") + ylab("TAC adjustment") + ggthemes::theme_few()  #+ ylab('TAC adjustment based on M. paradoxus biomass')```{r lam_gam}
ggsave(here("mods", "figs", "gamma.png"), width = 6, height = 5)
````

#figure([
#box(image("00-Namibian_hake_model_2024_files/figure-typst/fig-gamma-1.svg"))
], caption: figure.caption(
position: bottom, 
[
TAC adjustments for different levels of current biomass of M. paradoxus biomass.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-gamma>


```r
#---Plot the adjustment factors for different reactivities to M. paradoxus biomass----
df01 <- tibble()
df01 <- rbind(df01, tibble(relbiom = 1:100/100, lambda = rep(0.5, 100)), tibble(relbiom = 1:100/100,
    lambda = rep(1, 100)), tibble(relbiom = 1:100/100, lambda = rep(1.5, 100)), tibble(relbiom = 1:100/100,
    lambda = rep(2, 100))) |>
    mutate(adjustment = if_else(relbiom < 0.5, relbiom^lambda/0.5^lambda, 1), lambda = as.factor(lambda))
df01 |>
    ggplot(aes(y = adjustment, x = relbiom, color = lambda)) + geom_line() + xlab("Spawning biomass relative Bmsy") +
    ylab("TAC adjustment") + ggthemes::theme_few()  #+ ylab('TAC adjustment based on M. paradoxus biomass')
ggsave(here("mods", "figs", "lambda.png"), width = 6, height = 5)
```

#figure([
#box(image("00-Namibian_hake_model_2024_files/figure-typst/fig-lambda-1.svg"))
], caption: figure.caption(
position: bottom, 
[
TAC adjustments for different levels of current biomass relative to Bmsy
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-lambda>


= Notes
<notes>
Throughout the weeks, we scrutinized data inputs and found a couple of inconsistencies. For example, data were provided for surveys in 2019 yet there were no observations in that year. Similarly, the abundance-at-length data from the surveys failed to identify strong periods of persistence by cohorts.

Ignore 7-vessel CPUE data set.
