import { useAppStore } from "../store/appStore";
import type { Locale, Translations } from "./types";
import { en } from "./en";
import { nl } from "./nl";
import { de } from "./de";
import { es } from "./es";
import { fr } from "./fr";
import { pt } from "./pt";
import { zh } from "./zh";

export type { Locale, Translations };

const locales: Record<Locale, Translations> = { en, nl, de, es, fr, pt, zh };

/** Detect the best matching locale from the browser/OS language setting. */
export function detectSystemLocale(): Locale {
  const lang = navigator.language.toLowerCase();
  if (lang.startsWith("nl")) return "nl";
  if (lang.startsWith("de")) return "de";
  if (lang.startsWith("es")) return "es";
  if (lang.startsWith("fr")) return "fr";
  if (lang.startsWith("pt")) return "pt";
  if (lang.startsWith("zh")) return "zh";
  return "en";
}

/** React hook — returns a t() function bound to the current locale. Re-renders on locale change. */
export function useT(): (key: keyof Translations) => string {
  const locale = useAppStore((s) => s.locale);
  const dict = locales[locale] ?? en;
  return (key: keyof Translations) => dict[key];
}

/** Non-reactive helper for use outside components (e.g. Rust report command payload). */
export function getTranslations(locale: Locale): Translations {
  return locales[locale] ?? en;
}

export { en, nl, de, es, fr, pt, zh };
