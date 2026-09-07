# asariR command reference

本页放置README中不适合展开的参数级说明。普通用户先阅读
[`README.md`](../README.md)，只在需要调整算法参数时查询本页。

## `asari_process()`

```r
asari_process(
  input,
  output = "asari-output",
  project_name = "asari_project",
  mode = NULL,
  ppm = NULL,
  multicores = NULL,
  rt_align = NULL,
  database_mode = NULL,
  recursive = FALSE,
  parameters = list()
)
```

- `input`：一个明硤mzML文件、mzML路径向量或目录。
- `output`：项目输出根目录。
- `project_name`：本次运行的名称。
- `mode`：`"pos"`或`"neg"`。
- `ppm`：m/z容差。
- `multicores`：并行工作数。
- `rt_align`：多样本时是否执行保留时间对齐。
- `database_mode`：`"memory"`、`"ondisk"`或`"auto"`。
- `recursive`：输入是目录时是否搜索子目录。
- `parameters`：不常用参数的具名list。

## 高级参数

当前默认值以运行时对象为准：

```r
defaults <- asari_default_parameters()
View(as.data.frame(defaults))
```

常调整的字段：

| 字段 | 含义 |
|---|---|
| `mass_range` | 参与处理的m/z范围 |
| `min_timepoints` | 质量轨迹的最少数据点 |
| `min_intensity_threshold` | 输入数据点最低强度 |
| `min_peak_height` | 色谱峰最低高度 |
| `signal_noise_ratio` | 推荐表使用的信噪比标准 |
| `gaussian_shape` | 推荐表使用的峰形标准 |
| `peak_area` | `"sum"`、`"auc"`或`"gauss"` |
| `reference` | 多样本对齐的参考样本 |
| `max_retention_shift` | 允许的最大RT偏移 |
| `output_feature_table` | 特征表文件名主体 |

覆盖示例：

```r
result <- asari_process(
  input = "/absolute/path/to/data",
  output = "/absolute/path/to/results",
  parameters = list(
    min_timepoints = 5L,
    min_intensity_threshold = 2000,
    min_peak_height = 50000,
    peak_area = "auc"
  )
)
```

## 特征表固定列

| 列 | 含义 |
|---|---|
| `id_number` | 本次运行中的特征ID |
| `mz` | 特征m/z |
| `rtime` | 峰顶保留时间，秒 |
| `rtime_left_base` | 左峰界，秒 |
| `rtime_right_base` | 右峰界，秒 |
| `parent_masstrack_id` | 所属质量轨迹ID |
| `peak_area` | 跨样本汇总峰面积 |
| `cSelectivity` | 色谱选择性 |
| `goodness_fitting` | 峰形拟合质量 |
| `snr` | 信噪比 |
| `detection_counts` | 检出该特征的样本数 |

第12列开始是各样本的强度。

## 公开函数查询

R安装包后的帮助页是最新参数来源：

```r
help(package = "asariR")
?asari_process
?asari_read_results
?asari_annotation
?asari_ms2
?asari_dashboard
?asari_raw_converter
```
