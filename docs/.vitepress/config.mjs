import { defineConfig } from 'vitepress'

// https://vitepress.dev/reference/site-config
export default defineConfig({
  title: " ",
  // description: "A VitePress Site",
  themeConfig: {
    logo: 'https://vitepress.dev/vitepress-logo-mini.svg',
    // https://vitepress.dev/reference/default-theme-config
    nav: [
      { text: 'Research', link: 'https://github.com' },
      { text: 'CellAgent', link: 'https:/github.com' }
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
