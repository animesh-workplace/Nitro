<template>
    <div>
        <v-chart
            autoresize
            class="chart"
            ref="recombinant"
            :loading="loaded"
            :option="options"
            @finished="testfunction"
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
        loaded: true,
        loader_option: {
            color: '#c23531',
            lineWidth: 3,
            fontSize: 14,
            text: 'Loading',
            fontWeight: 500,
            fontFamily: 'Averta',
        },
        options: {
            grid: {left: '5%', right: '5%', top: '20%', bottom: '20%'},
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
                data: [
                    '271',
                    '1293',
                    '4486',
                    '11020',
                    '29158',
                    '29903',
                    '29904',
                    '29904',
                    '29904',
                    '29904',
                    '29904',
                    '29904',
                    '29904',
                    '29904',
                ],
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
                    barWidth: '95%',
                    itemStyle: {
                        borderRadius: [3, 3, 3, 3],
                    },
                    label: {
                        formatter: (params) => {
                            console.log(params.dataIndex)
                            return 'K'
                        },
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
    methods: {
        testfunction() {
            console.log('finished')
        }
    },
    mounted() {
        this.$nextTick(() => {
            // console.log(this.$refs.recombinant.resize())
            if(this.$refs.recombinant.chart) {
                this.loaded = false
            }
            // console.log(this.$refs.recombinant.resize({width:'auto', height:'auto'}))
            // this.$refs.recombinant.resize()
        })
    },
}
</script>

<style scoped>
.chart {
    height: 35px;
    width: 100%;
}
</style>
