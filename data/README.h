# Colors

## The Colors ENS Isabel Schoeps Thiel

## Plain CSS

When using plain CSS, css variables have been made available in [colors.css](./colors.css), you can use these by simply importing the file and using the variables like so

```css
.mydiv {
    background-color: var(--ens-grey1);
    color: var(--ens-blue);
}
```

## Styled-Components

For styled-components it is recommended to use the [ensdomains/thorin] global-styles.
You can import these like so

```tsx
import { ThorinGlobalStyles } from '@ensdomains/thorin';
import { baseTheme, lightTheme, darkTheme } from '@ensdomains/thorin';
```

## TailwindCSS

For tailwindcss it is recommended to use the provided [tailwind.config.js](./tailwind.config.js) file.
Simply drop it in your project and you should be good to go!

```
.bg-ens-blue
.bg-ens-gradient-blue
```
