<template>
    <div>
        <v-chart
            v-show="loaded"
            class="chart"
            ref="recombinant"
            :loading="false"
            :option="options"
            :loading-options="loader_option"
        />
    </div>
</template>

<script>
import { use } from 'echarts/core'
import { GridComponent, LegendComponent, TooltipComponent } from 'echarts/components'
import { BarChart } from 'echarts/charts'
import VChart, { THEME_KEY } from 'vue-echarts'
import { CanvasRenderer } from 'echarts/renderers'

use([BarChart, GridComponent, CanvasRenderer, LegendComponent, TooltipComponent])

export default {
    data: () => ({
        loaded: false,
        loader_option: {
            color: '#c23531',
            lineWidth: 3,
            fontSize: 14,
            text: 'Loading',
            fontWeight: 500,
            fontFamily: 'Averta',
        },
        options: {
            xAxis: {
                type: 'category',
                show: false,
                // show: true,
                // axisLabel: {
                //     rotate: '90',
                // },
                axisLine: {
                    show: false,
                },
                axisTick: {
                    show: false,
                },
                data: ['271', '1293', '4486', '11020', '29158', '29903', '29904','29904', '29904', '29904', '29904', '29904', '29904', '29904'],
                position: 'top',
            },
            yAxis: {
                show: false,
            },
            series: [
                {
                    data: [
                        100,
                        {
                            value: 100,
                            itemStyle: {
                                color: '#a90000',
                            },
                        },
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        {
                            value: 100,
                            itemStyle: {
                                color: '#a90000',
                            },
                        },
                        100,
                        100,
                    ],
                    type: 'bar',
                    barWidth: '98%',
                    itemStyle: {
                        borderRadius: [3, 3, 3, 3],
                    },
                    label: {
                        formatter: 'K',
                        show: true,
                    },
                },
            ],
        },
    }),
    computed: {},
    components: {
        VChart,
    },
    provide: {
        [THEME_KEY]: 'light',
    },
    methods: {},
    mounted() {
        this.$nextTick(() => {
            // console.log(this.$refs.recombinant.resize())
            this.loaded = true
            this.$refs.recombinant.resize()
        })
    },
}
</script>

<style scoped>
.chart {
    height: 100%;
    /*width: 100%;*/
}
</style>
