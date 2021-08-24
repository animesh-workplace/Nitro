import SuperTokens from 'supertokens-website'

export default function ({ $axios, redirect }) {
    SuperTokens.addAxiosInterceptors($axios)
}
