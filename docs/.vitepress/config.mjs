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
      { text: 'Research', link: '#' },
      { text: 'CellAgent', link: 'http://cell.agent4science.cn/start/' }
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
      { icon: 'github', link: '#' }
    ]
  },
  vite: {
    ssr: {
      noExternal: ['vuetify']
    }
  }
})
