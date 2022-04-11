const defaultTheme = require('tailwindcss/defaultTheme')

const fontFamily = defaultTheme.fontFamily
fontFamily['sans'] = ['Lexend Deca', 'system-ui']

module.exports = {
    content: [
        './pages/**/*.vue',
        './layouts/**/*.vue',
        './nuxt.config.{js,ts}',
        './plugins/**/*.{js,ts}',
        './node_modules/tv-*/**/*.vue',
        './components/**/*.{js,vue,ts}',
    ],
    theme: {
        fontFamily: fontFamily,
    },
    plugins: [require('flowbite/plugin')],
}
