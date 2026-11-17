import { withMermaid } from 'vitepress-plugin-mermaid'

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
      { text: 'Guide',        link: '/guide/' },
      { text: 'Tutorials',    link: '/tutorials/' },
      { text: 'Cookbook',     link: '/cookbook/' },
      { text: 'Verification', link: '/verification/' },
      { text: 'Regression',   link: '/regression/' },
      { text: 'API',          link: '/api/' },
      { text: 'GitHub',       link: 'https://github.com/szaghi/GIDEON' },
    ],

    sidebar: {
      '/guide/': [
        { text: 'Introduction', items: [
          { text: 'About',        link: '/guide/' },
          { text: 'Installation', link: '/guide/installation' },
          { text: 'Quick Start',  link: '/guide/quickstart' },
        ] },
        { text: 'The library', items: [
          { text: 'Conventions',  link: '/guide/conventions' },
          { text: 'Operators',    link: '/guide/operators' },
          { text: 'Precision',    link: '/guide/precision' },
          { text: 'API',          link: '/guide/api' },
          { text: 'The generator', link: '/guide/generator' },
          { text: 'Performance',  link: '/guide/performance' },
          { text: 'Sources',      link: '/guide/sources' },
        ] },
        { text: 'Project', items: [
          { text: 'Contributing', link: '/guide/contributing' },
          { text: 'Changelog',    link: '/guide/changelog' },
        ] },
      ],
      '/tutorials/': [
        { text: 'Tutorials', items: [
          { text: 'Overview',                         link: '/tutorials/' },
          { text: 'T1 A derivative on a periodic grid', link: '/tutorials/fd-derivative' },
          { text: 'T2 Finite volumes, conservatively', link: '/tutorials/fv-flux-difference' },
          { text: 'T3 Upwind faces and their mirror',  link: '/tutorials/upwind' },
          { text: 'T4 The order at run time',          link: '/tutorials/run-time-order' },
          { text: 'T5 Gradient, divergence, curl',     link: '/tutorials/operators-3d' },
          { text: 'T6 On the GPU',                     link: '/tutorials/gpu' },
        ] },
      ],
      '/cookbook/': [
        { text: 'Cookbook', items: [
          { text: 'Overview',                        link: '/cookbook/' },
          { text: 'Choose an operator and its order', link: '/cookbook/choose' },
          { text: 'Ghost cells for a window',         link: '/cookbook/ghosts' },
          { text: 'Replace ADAM\'s FDV calls',        link: '/cookbook/from-adam' },
          { text: 'Averages and point values',        link: '/cookbook/deconvolution' },
          { text: 'A run-time order, fast',           link: '/cookbook/run-time-order' },
          { text: 'Single and quadruple precision',   link: '/cookbook/precision' },
          { text: 'Add an operator',                  link: '/cookbook/new-operator' },
        ] },
      ],
      '/verification/': [
        { text: 'Verification', items: [
          { text: 'Catalog', link: '/verification/' },
        ] },
        { text: 'Contracts', collapsed: false, items: [
          { text: 'C1 Rows by their meaning',   link: '/verification/rows' },
          { text: 'C2 Exactness and leading error', link: '/verification/exactness' },
          { text: 'C3 Mirrors and symmetries',  link: '/verification/properties' },
          { text: 'C4 Run-time dispatch',       link: '/verification/dispatch' },
          { text: 'C5 Equivalence with ADAM',   link: '/verification/adam' },
          { text: 'C6 Host and device',         link: '/verification/device' },
        ] },
        { text: 'Accuracy', collapsed: false, items: [
          { text: 'A1 Orders of accuracy',      link: '/verification/orders' },
          { text: 'A2 3-D operators',           link: '/verification/operators-3d' },
          { text: 'A3 Spectral properties',     link: '/verification/spectral' },
        ] },
      ],
      '/regression/': [
        { text: 'Regression', items: [
          { text: 'Overview',         link: '/regression/' },
          { text: 'Golden values',    link: '/regression/golden' },
          { text: 'Mutation testing', link: '/regression/mutation' },
          { text: 'Continuous integration', link: '/regression/ci' },
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
