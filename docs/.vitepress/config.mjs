import { defineConfig } from 'vitepress'

// https://vitepress.dev/reference/site-config
export default defineConfig({
  title: "CellAgent",
  // description: "A VitePress Site",
  themeConfig: {
    // https://vitepress.dev/reference/default-theme-config
    siteTitle: '',
    logo: '/logo.png',
    nav: [
      { text: 'Research', link: 'https://www.biorxiv.org/content/10.1101/2024.05.13.593861v1' },
      { text: 'CellAgent', link: 'http://cell.agent4science.cn/' }
    ],

    sidebar: [
      {
        text: 'Examples',
        items: [
          { text: 'Markdown Examples', link: '/markdown-examples' },
          { text: 'Runtime API Examples', link: '/api-examples' }
        ]
      }
    ],

    socialLinks: [
      { icon: 'github', link: 'https://github.com/lichaozhy/cell-agent-home' }
    ]
  },
  vite: {
    ssr: {
      noExternal: ['vuetify']
    }
  }
})
