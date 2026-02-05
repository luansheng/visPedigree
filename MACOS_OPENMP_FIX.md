# macOS R OpenMP 环境修复备忘录

**日期**: 2026-02-05 **问题**: 在 macOS 上，R 默认自带的 OpenMP动态库
(`libomp.dylib`) 是残缺的，缺少某些符号（如
`___kmpc_dispatch_deinit`），导致在使用 RcppArmadillo
进行多线程开发时出现 `symbol not found` 崩溃。

## 执行的操作

我们执行了 “外科手术式” 的替换，用 Homebrew 提供的完整版 `libomp` 替换了
R 自带的版本。

### 1. 备份原文件

``` bash
sudo mv /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libomp.dylib \
        /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libomp.dylib.bak
```

*如果未来需要恢复，只需把这个 `.bak` 文件改回原名即可。*

### 2. 创建链接

``` bash
sudo ln -s /opt/homebrew/opt/libomp/lib/libomp.dylib \
           /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libomp.dylib
```

这使得 R 在加载 OpenMP 时，实际上加载的是完整功能的 Homebrew 版本。

## 潜在副作用评估

**风险等级：极低**

1.  **对 R 的其他包有影响吗？**
    - 绝大多数情况下是 **正面的**。其他依赖 OpenMP 的 R 包（如
      `data.table` 等）也会因此受益，获得更稳定的多线程支持。
    - Homebrew 的 `libomp` 是 LLVM 项目的标准实现，向下兼容性很好。
2.  **R 升级后需要重做吗？**
    - **是的**。如果您安装了 R 的新版本（例如从 4.5.2 升级到 4.6.0），R
      安装程序会写入新的
      `libomp.dylib`，覆盖掉我们的链接。到时候您可能再次遇到同样的报错，只需**重新执行上述两步**即可。
3.  **对系统其他软件有影响吗？**
    - **没有**。我们只修改了 R
      框架内部的文件，不影响系统自带的库或其他软件。

## 恢复方法 (如果有任何异常)

如果发现 R 崩溃或无法启动，运行以下命令恢复原状：

``` bash
sudo rm /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libomp.dylib
sudo mv /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libomp.dylib.bak \
        /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libomp.dylib
```
