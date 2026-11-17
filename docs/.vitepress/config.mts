import { existsSync, readFileSync } from 'node:fs'
import { resolve } from 'node:path'
import { withMermaid } from 'vitepress-plugin-mermaid'

// formal's sidebar of the API pages (docs/api/_sidebar.json), when the API has been generated
const apiSidebarFile = resolve(__dirname, '../api/_sidebar.json')
const apiSidebar = existsSync(apiSidebarFile) ? JSON.parse(readFileSync(apiSidebarFile, 'utf-8')) : []

export default withMermaid({
  title: 'GIDEON',
  description: 'Generated Interpolations and Derivatives, Exact to Order N: finite-difference and finite-volume stencil operators for HPC',
  base: '/GIDEON/',
  srcExclude: ['ford.md', '**/_tables/**', 'guide/_operators-table.md', 'examples/**'],

  head: [
    ['meta', { name: 'theme-color', content: '#f3ede1', media: '(prefers-color-scheme: light)' }],
    ['meta', { name: 'theme-color', content: '#0e1116', media: '(prefers-color-scheme: dark)' }],
    // the theme-aware videos of the pages: <video class="light-only"> / <video class="dark-only">
    ['style', {}, '.dark-only{display:none}.dark .dark-only{display:block}.dark .light-only{display:none}'],
  ],

  markdown: {
    math: true,
    languages: ['fortran-free-form', 'fortran-fixed-form'],
    languageAlias: {
      fortran: 'fortran-free-form',
      f90: 'fortran-free-form',
      f03: 'fortran-free-form',
      f08: 'fortran-free-form',
    },
  },

  themeConfig: {
    siteTitle: 'GIDEON',
    nav: [
      { text: 'Start',        link: '/guide/' },
      { text: 'Tutorials',    link: '/tutorials/' },
      { text: 'Cookbook',     link: '/cookbook/' },
      { text: 'Verification', link: '/verification/' },
      { text: 'Regression',   link: '/regression/' },
      { text: 'API',          link: '/api/' },
      { text: 'GitHub',       link: 'https://github.com/szaghi/GIDEON' },
    ],

    // One sidebar for the whole site, in reading order: it is also the chain of the Previous / Next links at the foot
    // of every page, so reading straight through goes start -> tutorials -> recipes -> reference -> verification.
    // A collapsed group opens by itself on its own pages. scripts/check_docs.py fails when a page is missing here.
    sidebar: {
      // the module pages written by formal (docs:build runs it first), with its own sidebar
      '/api/': [
        { text: 'Back to the documentation', link: '/guide/' },
        { text: 'API', items: [{ text: 'Overview', link: '/api/' }] },
        ...apiSidebar,
      ],
      '/': [
        { text: 'Start here', items: [
          { text: 'About',        link: '/guide/' },
          { text: 'Installation', link: '/guide/installation' },
          { text: 'Quick Start',  link: '/guide/quickstart' },
          { text: 'Conventions',  link: '/guide/conventions' },
        ] },
        { text: 'Tutorials', collapsed: false, items: [
          { text: 'All tutorials',                       link: '/tutorials/' },
          { text: 'T1 A derivative on a periodic grid',  link: '/tutorials/fd-derivative' },
          { text: 'T2 Finite volumes, conservatively',   link: '/tutorials/fv-flux-difference' },
          { text: 'T3 Upwind faces and their mirror',    link: '/tutorials/upwind' },
          { text: 'T4 The order at run time',            link: '/tutorials/run-time-order' },
          { text: 'T5 Gradient, divergence, curl',       link: '/tutorials/operators-3d' },
          { text: 'T6 On the GPU',                       link: '/tutorials/gpu' },
        ] },
        { text: 'Cookbook', collapsed: false, items: [
          { text: 'All recipes',                      link: '/cookbook/' },
          { text: 'Choose an operator and its order', link: '/cookbook/choose' },
          { text: 'Ghost cells for a window',         link: '/cookbook/ghosts' },
          { text: 'Averages and point values',        link: '/cookbook/deconvolution' },
          { text: 'A run-time order, fast',           link: '/cookbook/run-time-order' },
          { text: 'Single and quadruple precision',   link: '/cookbook/precision' },
          { text: 'Replace ADAM\'s FDV calls',        link: '/cookbook/from-adam' },
          { text: 'Add an operator',                  link: '/cookbook/new-operator' },
        ] },
        { text: 'Reference', collapsed: false, items: [
          { text: 'Operators',     link: '/guide/operators' },
          { text: 'API',           link: '/guide/api' },
          { text: 'Precision',     link: '/guide/precision' },
          { text: 'The generator', link: '/guide/generator' },
          { text: 'Performance',   link: '/guide/performance' },
        ] },
        { text: 'Verification', collapsed: true, items: [
          { text: 'Verification catalog', link: '/verification/' },
          { text: 'C1 Rows by their meaning',       link: '/verification/rows' },
          { text: 'C2 Exactness and leading error', link: '/verification/exactness' },
          { text: 'C3 Mirrors and symmetries',      link: '/verification/properties' },
          { text: 'C4 Run-time dispatch',           link: '/verification/dispatch' },
          { text: 'C5 Equivalence with ADAM',       link: '/verification/adam' },
          { text: 'C6 Host and device',             link: '/verification/device' },
          { text: 'A1 Orders of accuracy',          link: '/verification/orders' },
          { text: 'A2 3-D operators',               link: '/verification/operators-3d' },
          { text: 'A3 Spectral properties',         link: '/verification/spectral' },
        ] },
        { text: 'Regression', collapsed: true, items: [
          { text: 'Regression overview',    link: '/regression/' },
          { text: 'Golden values',          link: '/regression/golden' },
          { text: 'Mutation testing',       link: '/regression/mutation' },
          { text: 'Continuous integration', link: '/regression/ci' },
        ] },
        { text: 'Project', collapsed: true, items: [
          { text: 'Sources',      link: '/guide/sources' },
          { text: 'Contributing', link: '/guide/contributing' },
          { text: 'Changelog',    link: '/guide/changelog' },
        ] },
      ],
    },

    socialLinks: [
      { icon: 'github', link: 'https://github.com/szaghi/GIDEON' },
    ],

    search: {
      provider: 'local',
    },

    footer: {
      message: 'Released under your choice of GPL-3.0-or-later, BSD-2-Clause, BSD-3-Clause or MIT.',
      copyright: 'Copyright © 2026 Stefano Zaghi',
    },
  },

  mermaid: {},

  vite: {
    // es2022: mermaid 11.16+ emits syntax that esbuild cannot lower to vite's default es2020 target
    build: {
      target: 'es2022',
    },
  },
})
