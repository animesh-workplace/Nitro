<template>
    <div>
        <!-- <Tutorial /> -->
        <!--         <div class="grid grid-cols-3 gap-4">
            <ChartsRecombinant class="col-start-2 col-span-2" />
        </div> -->
        <RtableTable />
        <!-- <nuxt-link to="/test">LINK</nuxt-link> -->
    </div>
</template>

<script>
import SuperTokens from 'supertokens-website'

export default {
    name: 'IndexPage',
    mounted() {
        this.$nextTick(async () => {
            // SuperTokens.init({
            //     apiDomain: 'http://localhost:3000',
            //     apiBasePath: '/auth',
            // })
            // console.log(await SuperTokens.doesSessionExist())
            // console.log('getting user information')
            // console.log(await SuperTokens.getUserId())
            console.log(this.$axios)
        })
    },
}
</script>
