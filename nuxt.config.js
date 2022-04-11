export default {
    css: [],
    plugins: [],
    components: true,
    axios: { baseURL: '/' },
    modules: ['@nuxtjs/axios'],
    build: {
        transpile: [/echarts/, /zrender/],
    },
    buildModules: ['@nuxtjs/tailwindcss', '@nuxtjs/composition-api/module', '@nuxtjs/google-fonts'],
    head: {
        title: 'Recombinant Finder',
        htmlAttrs: {
            lang: 'en',
        },
        meta: [
            { charset: 'utf-8' },
            { name: 'viewport', content: 'width=device-width, initial-scale=1' },
            { hid: 'description', name: 'description', content: '' },
            { name: 'format-detection', content: 'telephone=no' },
        ],
        link: [{ rel: 'icon', type: 'image/x-icon', href: '/favicon.ico' }],
    },
    googleFonts: {
        preload: true,
        prefetch: true,
        download: true,
        display: 'swap',
        preconnect: true,
        overwriting: false,
        families: {
            'Lexend+Deca': ['100', '200', '300', '400', '500', '600', '700', '800', '900'],
        },
    },
}
