<template>
  <v-sheet class="my-16 position-relative" border="sm" rounded>
    <v-avatar
      size="33"
      class="chat-bubble-subject position-absolute border"
    >
      <slot name="subject">
        <v-icon :icon="mdiAccount"></v-icon>
      </slot>
    </v-avatar>

    <div class="pa-2">
      <slot>
        <div>content</div>
      </slot>
    </div>
  </v-sheet>

</template>

<script setup>
import { mdiAccount } from '@mdi/js'
</script>

<style>
.chat-bubble-subject {
  top: -40px;
}
</style>
